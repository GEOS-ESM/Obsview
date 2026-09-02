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
from .stats.aggregate import str_to_datetime, aggregate_stats, aggregate_pass_binned, aggregate_fail_binned
from .plotting.statsplot import plot_stats
from .plotting.spatialcoverage import plot_coverage
from .plotting.timeseries import plot_series
from .plotting.radmon import plot_radmon
from. plotting.compare_statsplot import plot_compare_stats




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

#TODO: Make module for searching experiment directories and finding specific tar files for comparing experiments
def experiment_tar_dir(base_path: str, expid: str, dt: str) -> str:
    """
    Build the tarball directory for one experiment and a given datetime.

    Layout: <base_path>/<expid>/jedi/obs/Y<YYYY>/M<MM>/

    Parameters
    ----------
    base_path : str
        Absolute path to the directory containing experiment folders.
    expid : str
        Experiment id, e.g. 'j54rp1'.
    dt : str
        Expected string format: 'YYYYMMDDHH', example: '2026010100'
    """
    dt = str_to_datetime(dt)
    
    if not os.path.isabs(base_path):
        raise ValueError(f"base_path must be absolute: {base_path!r}")

    year_dir = f"Y{dt.year:04d}"
    month_dir = f"M{dt.month:02d}"
    tar_dir = os.path.join(base_path, expid, "jedi", "obs", year_dir, month_dir)

    if not os.path.isdir(tar_dir):
        raise FileNotFoundError(f"Experiment tar directory not found: {tar_dir!r}")
    return tar_dir



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
    filename = "python/data/ODS files/radiance/j54rp1.diag_atms_n20.20260101_00z.ods"
    tarpath = "/Users/ltrayano/Desktop/Obsview/Obsview/python/data/IODA files"
    instrument = "atms_n20"
    varname = "brightnessTemperature"
    kx = 920
    starttime = "2026010100"
    endtime = "2026013118"

    reader = ODSReader()
    data = reader.read(filename, varname, kx)
    
    #masking
    val_mask = fill_val_mask(data)

    #ts_plot = plot_series(ts)
    #filtering
    valid_data = apply_filter(data, val_mask)

    

    
    #QC masking
    pass_mask = qc_pass_mask(valid_data)
    fail_mask = qc_fail_mask(valid_data)


    pass_data = apply_filter(valid_data, pass_mask)
    

    #calculate job, joa, esigo, esigb
    pass_data = calc_derived(pass_data)                 #Only calculate variables for QC = 0 data

    #binning
    pass_data_binned = create_bins(data)

    radmon_plot = plot_radmon(pass_data_binned)
    plt.show()

    

    





    # filepath1 = "python/data/IODA files/j54rp1/jedi/obs/Y2026/M01/j54rp1.jedi_hofx.20260101_03z/atms_n20.20260101T030000Z.nc4"
    # filepath2 = "python/data/IODA files/j54rp2/jedi/obs/Y2026/M01/j54rp2.jedi_hofx.20260101_03z/atms_n20.20260101T030000Z.nc4"
    # instrument = "atms_n20"
    # varname = "brightnessTemperature"
    # kx = 920
    # starttime = "2026010100"
    # endtime = "2026013118"

    # exp1 = "j54rp1"
    # exp2 = "j54rp2"

    # reader = IODAReader()
    # #File 1
    # data = reader.read(filepath1, varname, kx)
    # val_mask = fill_val_mask(data)    
    #             #filtering
    # valid_data = apply_filter(data, val_mask)         
    #             #QC masking
    # pass_mask = qc_pass_mask(valid_data)       
    # pass_data = apply_filter(valid_data, pass_mask)

    #             #calculate job, joa, esigo, esigb
    # pass_data = calc_derived(pass_data)                 
    #             #binning
    # pass_data1_binned = create_bins(pass_data)          
    #             #stats
    # pass_stats1_binned = calculate_stats(pass_data1_binned)

    # #File 2
    # data = reader.read(filepath2, varname, kx)
    # val_mask = fill_val_mask(data)    
    #             #filtering
    # valid_data = apply_filter(data, val_mask)         
    #             #QC masking
    # pass_mask = qc_pass_mask(valid_data)       
    # pass_data = apply_filter(valid_data, pass_mask)

    #             #calculate job, joa, esigo, esigb
    # pass_data = calc_derived(pass_data)                 
    #             #binning
    # pass_data2_binned = create_bins(pass_data)          
    #             #stats
    # pass_stats2_binned = calculate_stats(pass_data2_binned)

    # compare_plot = plot_compare_stats(pass_data1_binned, pass_stats1_binned, pass_stats2_binned)
    # plt.show()











    

    
    




if __name__ == '__main__':
    main()