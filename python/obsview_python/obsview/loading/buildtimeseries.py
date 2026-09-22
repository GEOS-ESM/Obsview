#Module containing functions to build a time series object from ODS, IODA, and TAR files
import os
import glob
import tarfile
import tempfile
from dataclasses import replace
from typing import List, Optional
from pathlib import Path
from functools import partial
from concurrent.futures import ProcessPoolExecutor

from ..loading.timeseriesdata import TimeSeriesData
from ..loading.iodareader import IODAReader
from ..loading.odsreader import ODSReader
from ..loading.tarreader import find_instrument_member, extract_member_path, _process_single_tar
from ..processing.masking import fill_val_mask, qc_pass_mask, qc_fail_mask
from ..processing.filtering import apply_filter
from ..processing.derived import calc_derived
from .. processing.binning import create_bins
from ..stats.calc_stats import calculate_stats
from ..stats.aggregate import str_to_datetime




#Function to loop over IODA files and create a timeseries data object
def build_ts_from_ioda(filenames: List[str], varname: str, kx: int) -> TimeSeriesData:
    reader = IODAReader()
    records = []        #List to append relevant contents of each file to

    for fname in filenames:
        data = reader.read(fname, varname, kx)
        #masking
        val_mask = fill_val_mask(data)
    
        #filtering
        valid_data = apply_filter(data, val_mask)
        
        #QC masking
        pass_mask = qc_pass_mask(valid_data)
        fail_mask = qc_fail_mask(valid_data)
        pass_data = apply_filter(valid_data, pass_mask)
        fail_data = apply_filter(valid_data, fail_mask)

        #calculate job, joa, esigo, esigb
        pass_data = calc_derived(pass_data)                 #Only calculate variables for QC = 0 data

        #binning
        pass_data_binned = create_bins(pass_data)
        fail_data_binned = create_bins(fail_data)
        #stats
        pass_stats_binned = calculate_stats(pass_data_binned)      #Only calculate stats on QC = 0 data

        #append objects to records list
        records.append(
            (data.datetime, pass_stats_binned, pass_data_binned, fail_data_binned)
        )

    records.sort(key=lambda r: r[0])        #Sort chronologically

    #Append datetimes
    datetimes  = [r[0] for r in records]
    pass_stats = [r[1] for r in records]
    pass_data  = [r[2] for r in records]
    fail_data  = [r[3] for r in records]

    obj = TimeSeriesData(
        datetimes = datetimes,
        pass_stats = pass_stats,
        pass_data = pass_data,
        fail_data = fail_data 
    )
    return obj
    ...

def build_ts_from_ods(filenames: List[str], varname, kx, start_time: str, end_time: str) -> TimeSeriesData:
    reader = ODSReader()
    records = []        #List to append relevant contents of each file to
    starttime = str_to_datetime(start_time)
    endtime = str_to_datetime(end_time)

    for fname in filenames:
        data = reader.read(fname, varname, kx)
        #masking
        val_mask = fill_val_mask(data)
    
        #filtering
        valid_data = apply_filter(data, val_mask)
        
        #QC masking
        pass_mask = qc_pass_mask(valid_data)
        fail_mask = qc_fail_mask(valid_data)
        pass_data = apply_filter(valid_data, pass_mask)
        fail_data = apply_filter(valid_data, fail_mask)

        #calculate job, joa, esigo, esigb
        pass_data = calc_derived(pass_data)                 #Only calculate variables for QC = 0 data

        #binning
        pass_data_binned = create_bins(pass_data)
        pass_data_binned = replace(pass_data_binned, ts_range = [starttime, endtime])
        fail_data_binned = create_bins(fail_data)
        #stats
        pass_stats_binned = calculate_stats(pass_data_binned)      #Only calculate stats on QC = 0 data

        #append objects to records list
        records.append(
            (data.datetime, pass_stats_binned, pass_data_binned, fail_data_binned)
        )

    records.sort(key=lambda r: r[0])        #Sort chronologically

    #Append datetimes
    datetimes  = [r[0] for r in records]
    pass_stats = [r[1] for r in records]
    pass_data  = [r[2] for r in records]
    fail_data  = [r[3] for r in records]

    obj = TimeSeriesData(
        datetimes = datetimes,
        pass_stats = pass_stats,
        pass_data = pass_data,
        fail_data = fail_data 
    )
    return obj
    ...

def build_ts_from_tar(tar_path: str, instrument: str, varname: str, kx: int, start_time: str, end_time: str) -> TimeSeriesData:
    tar_file_paths = sorted(glob.glob(os.path.join(tar_path, "*.tar")))
    starttime = str_to_datetime(start_time)
    endtime = str_to_datetime(end_time)

    # Freeze constant parameters for pool mapping
    worker_fn = partial(_process_single_tar,instrument=instrument,varname=varname,kx=kx,starttime=starttime,endtime=endtime)

    # Execute tasks across worker processes
    with ProcessPoolExecutor() as executor:
        # executor.map preserves the original sorted order of tar_file_paths
        results = executor.map(worker_fn, tar_file_paths)

    # Filter out None results (tar archives missing the instrument member)
    records = [r for r in results if r is not None]
    datetimes  = [r[0] for r in records]
    pass_stats = [r[1] for r in records]
    pass_data  = [r[2] for r in records]
    fail_data  = [r[3] for r in records]

    obj = TimeSeriesData(
        datetimes = datetimes,
        pass_stats = pass_stats,
        pass_data = pass_data,
        fail_data = fail_data 
    )
    return obj
    