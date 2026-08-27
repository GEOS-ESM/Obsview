#Main python script for now, running without command line arguments(to be changed later)
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional
import glob
import os
import tarfile
import tempfile

from .loading.odsreader import ODSReader
from .loading.iodareader import IODAReader
from .loading.timeseriesdata import TimeSeriesData
from .processing.masking import fill_val_mask, qc_pass_mask, qc_fail_mask
from .processing.filtering import apply_filter
from .processing.derived import calc_derived
from .processing.binning import create_bins
from .stats.calc_stats import calculate_stats
from .stats.aggregate import aggregate_stats, aggregate_pass_binned, aggregate_fail_binned
from .plotting.statsplot import plot_stats
from .plotting.spatialcoverage import plot_coverage
from .plotting.timeseries import plot_series
from .plotting.radmon import plot_radmon




#Function to find a specific .nc4 file inside a tarball of a specified instrument(e.g. atms_n20)
def find_instrument_member(tar: tarfile.TarFile, instrument: str) -> Optional[tarfile.TarInfo]:
    members = tar.getmembers()
    #Loop over all members of a tarball
    for member in members:
        base = os.path.basename(member.name)
        
        if base.startswith(instrument) and base.endswith("nc4"):
            return member
         
    return None
    ...

def extract_member_path(tar: tarfile.TarFile, member: tarfile.TarInfo, dest_dir: str) -> str:
     # Reject non-regular files (symlinks, hardlinks, devices, dirs-as-files).
    if not member.isfile():
        raise ValueError(f"Refusing to extract non-regular member: {member.name!r}")

    member_path = os.path.join(dest_dir, member.name)

    # Extract just this member. On Python 3.12+, filter='data' adds another
    # layer of protection; guard so it still works on older versions.
    try:
        tar.extract(member, path=dest_dir, filter="data")  # py3.12+
    except TypeError:
        tar.extract(member, path=dest_dir)                 # older Python

    extracted = os.path.abspath(member_path)
    return extracted

    ...



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


def build_ts_from_tar(tar_path: str, instrument: str, varname: str, kx: int) -> TimeSeriesData:
    tar_file_paths = sorted(glob.glob(os.path.join(tar_path, "*.tar")))
    records = []
    reader = IODAReader()

    #Loop through tar files
    for tar_file_path in tar_file_paths:
        with tarfile.open(tar_file_path, mode = "r:*") as tar:
            member = find_instrument_member(tar, instrument)
            if member == None:      #If instrument file is missing...
                continue
            with tempfile.TemporaryDirectory() as tmp:
                nc_path = extract_member_path(tar, member, tmp)
                data = reader.read(nc_path,varname, kx)
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


def build_ts_from_ods(filenames: List[str], varname, kx) -> TimeSeriesData:
    reader = ODSReader()
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


def make_map_plot(filename: str, varname: str, kx: int) -> None:
    #IODA file
    reader = IODAReader()
    data = reader.read(filename, varname, kx)
    

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

    #plotting
    coverage_plot = plot_coverage(pass_data_binned,fail_data_binned,14)
    plt.show()
    ...



def main() -> None:
    filename = "python/data/IODA files/j54rp1.jedi_hofx.20260101_03z/atms_n20.20260101T030000Z.nc4"
    tarpath = "/Users/ltrayano/Desktop/Obsview/Obsview/python/data/IODA files"
    instrument = "atms_n20"
    varname = "brightnessTemperature"
    kx = 920
    starttime = "2026010100"
    endtime = "2026013118"

    reader = IODAReader()
    data = reader.read(filename, varname, kx)
    

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

    radmon_plot = plot_radmon(pass_data_binned)
    plt.show()












    # ts = build_ts_from_tar(tarpath, instrument, varname, kx)
    
    #ts_plot = plot_series(ts)
    
    # ag_stats = aggregate_stats(ts,starttime, endtime)
    # ag_pass_binned = aggregate_pass_binned(ts, starttime, endtime)
    # ag_fail_binned = aggregate_fail_binned(ts, starttime, endtime)

    # stats_plot = plot_stats(ag_pass_binned,ag_fail_binned,ag_stats)
    # plt.show()


    
    




if __name__ == '__main__':
    main()