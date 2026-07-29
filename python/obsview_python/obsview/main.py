#Main python script for now, running without command line arguments(to be changed later)
import numpy as np
import matplotlib.pyplot as plt
from typing import List
import glob

from .loading.odsreader import ODSReader
from .loading.iodareader import IODAReader
from .loading.timeseriesdata import TimeSeriesData
from .processing.masking import valid_mask, fill_val_mask, valid_latlon_mask, qc_pass_mask, qc_fail_mask
from .processing.filtering import apply_filter
from .processing.derived import calc_derived
from .processing.binning import create_channel_bins
from .stats.calc_stats import calculate_stats
from .plotting.statsplot import plot_stats
from .plotting.spatialcoverage import plot_coverage
from .plotting.timeseries import plot_series





#Function to loop over IODA files and create a timeseries data object
def build_time_series(filenames: List[str]) -> TimeSeriesData:
    reader = IODAReader()
    records = []        #List to append relevant contents of each file to

    for fname in filenames:
        data = reader.read(fname)
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
        pass_data_binned = create_channel_bins(pass_data)
        fail_data_binned = create_channel_bins(fail_data)
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


def make_map_plot(filename: str) -> None:
    #IODA file
    reader = IODAReader()
    data = reader.read(filename)
    

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
    pass_data_binned = create_channel_bins(pass_data)
    fail_data_binned = create_channel_bins(fail_data)

    #plotting
    coverage_plot = plot_coverage(pass_data_binned,fail_data_binned,14)
    plt.show()
    ...



def main() -> None:

   filenames = sorted(glob.glob("python/data/IODA files/nc4/atms_n20.*.nc4"))

    
   ts = build_time_series(filenames)
   series_plot = plot_series(ts)
   plt.show()
   ...

    


if __name__ == '__main__':
    main()