#Main python script for now, running without command line arguments(to be changed later)
import time
import glob
import yaml
from pathlib import Path
import matplotlib.pyplot as plt
from dataclasses import replace

from .loading.buildtimeseries import build_ts_from_ods, build_ts_from_tar
from .stats.aggregate import str_to_datetime, aggregate_stats, aggregate_pass_binned, aggregate_fail_binned
from .plotting.statsplot import plot_stats
from .plotting.spatialcoverage import plot_coverage
from .plotting.timeseries import plot_series
from .plotting.radmon import plot_radmon
from. plotting.compare_statsplot import plot_compare_stats
from .plotting.purple import plot_purple









def main() -> None:
    config_path = Path(__file__).parent / 'config.yml'
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    

    ctl = config['ctl']
    ctl_path = config['ctl_path']
    exp_path = config['exp_path']
    instrument = config['instrument']
    varname = config['varname']
    kx = config['kx']
    starttime = str(config['starttime'])
    endtime = str(config['endtime'])

    filenames_ctl = sorted(glob.glob(f"python/data/ODS files/{ctl}/{ctl}.diag_atms_n20.*.ods"))

    print("Loading ODS...")
    start_time = time.perf_counter()
    ts_ctl = build_ts_from_ods(filenames_ctl,varname,kx, starttime, endtime)    
    end_time = time.perf_counter()
    print("ODS files loaded")
    print(f"Task time: {end_time - start_time} seconds")
    
    
    print("Loading IODA...")
    start_time = time.perf_counter()
    ts_exp = build_ts_from_tar(exp_path,instrument,varname,kx,starttime, endtime)
    end_time = time.perf_counter()
    print("IODA files loaded")
    print(f"Task time: {end_time - start_time} seconds")

    ctl_ag_bins = aggregate_pass_binned(ts_ctl,starttime,endtime)
    ctl_ag_stats = aggregate_stats(ts_ctl, starttime, endtime)
    exp_ag_stats = aggregate_stats(ts_exp, starttime, endtime)
    

    #stats_plot = plot_stats(ag_pass, ag_fail, ag_stats)

    #series_plot = plot_series(ts_exp, channel=9)

    purple_plot = plot_purple(ctl_ag_bins,ctl_ag_stats,exp_ag_stats, "std_omb")
    plt.show()
    




if __name__ == '__main__':
    main()