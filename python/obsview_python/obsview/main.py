#Main python script for now, running without command line arguments(to be changed later)
import matplotlib.pyplot as plt

from .loading.iodareader import IODAReader
from .processing.masking import valid_mask, fill_val_mask, qc_pass_mask, qc_fail_mask
from .processing.filtering import apply_filter
from .processing.derived import calc_derived
from .processing.binning import create_channel_bins
from .stats.calc_stats import calculate_stats
from .plotting.statsplot import plot_stats


filename = "python/amsua_metop-b.20260125T150000Z.nc4"
#@timer
def main() -> None:
    
    #IODA file
    reader = IODAReader()
    data = reader.read(filename)

    #masking
    valid_mask = fill_val_mask(data)
    
    #filtering
    valid_data = apply_filter(data, valid_mask)

    #QC masking
    pass_mask = qc_pass_mask(valid_data)
    fail_mask = qc_fail_mask(valid_data)


    pass_data = apply_filter(valid_data, pass_mask)
    fail_data = apply_filter(valid_data, fail_mask)

    #calculate job, joa, esigo, esigb
    pass_data = calc_derived(pass_data)

    #binning
    pass_data_binned = create_channel_bins(pass_data)
    fail_data_binned = create_channel_bins(fail_data)
    #stats
    data_stats = calculate_stats(pass_data_binned)

    #plotting
    panel_plot = plot_stats(pass_data_binned,fail_data_binned, data_stats)
    plt.show()




if __name__ == '__main__':
    main()