#Running script to execute code based on config settings
"""
Orchestration layer. Takes a validated Config and decides WHAT to run:
- resolves file locations from directory templates + time range,
- loops over experiments and targets,
- dispatches to the enabled plots,
- honors runtime settings (verbose, parallel) and output settings.

It calls into the existing pipeline/plotting modules but contains no
pipeline math itself.
"""

###########
#TODO: Reconcile start and endtime helper functions(i.e. have the datetime conversion happen once)

import os
import glob
import logging
from datetime import datetime, timedelta, timezone
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List
from dataclasses import replace

import matplotlib.pyplot as plt

from .config import Config, Experiment, Target

# Time series construction + plotting (existing modules).
from .utils import _parse_time
from .loading.timeseriesdata import TimeSeriesData
from .loading.search import discover_ods_files, discover_ioda_files
from .loading.buildtimeseries import build_ts_from_ods, build_ts_from_tar
from .processing.binning import BinnedData
from .stats.statisticsdata import StatisticsData
from .stats.aggregate import aggregate_stats, aggregate_pass_binned, aggregate_fail_binned
from .plotting.statsplot import plot_stats
from .plotting.timeseries import plot_series
from .plotting.spatialcoverage import plot_coverage
from .plotting.purple import plot_purple


log = logging.getLogger("obsview")


# ============================================================
# Entry point
# ============================================================

def run(cfg: Config) -> None:
    #_configure_logging(cfg)

    #log.info("Starting run: file_type=%s, %d experiment(s), %d target(s)",
             #cfg.file_type, len(cfg.experiments), len(cfg.targets))

    start = _parse_time(cfg.selection.time_range.start)
    end = _parse_time(cfg.selection.time_range.end)

    # Outer loop: each instrument/obtype the user asked to analyze.
    for target in cfg.selection.targets:
        
        #log.info("Target: instrument=%s obtype=%s varname=%s (%s levels)",
                 #target.instrument, target.obtype, varname, lev_type)

        # Load each experiment's time series for this target.
        ts_by_exp: Dict[str, TimeSeriesData] = _load_all_experiments(cfg, target, start, end)
        ag_by_exp = _aggregate_experiments(cfg,ts_by_exp, start, end)


        # Dispatch to enabled plots
        _make_enabled_plots(cfg, target, ts_by_exp, ag_by_exp, start, end)
        plt.show()

    #log.info("Run complete.")

# ============================================================
# Runtime settings
# ============================================================

# def _configure_logging(cfg: Config) -> None:
#     """Verbose flag controls log level; optional log file."""
#     level = logging.DEBUG if cfg.run.verbose else logging.INFO
#     handlers = [logging.StreamHandler()]
#     if getattr(cfg.run, "log_file", None):
#         handlers.append(logging.FileHandler(cfg.run.log_file))
#     logging.basicConfig(
#         level=level,
#         format="%(asctime)s [%(levelname)s] %(message)s",
#         handlers=handlers,
#     )


# ============================================================
# File location: templates + time range -> concrete directories
# ============================================================

def _resolve_tar_dir(cfg: Config, exp: Experiment, dt: datetime) -> str:
    """
    Fill the directory template from config with experiment + time values.

    Template example (from YAML):
      "{base}/{expid}/jedi/obs/Y{year}/M{month}"
    """
    template = cfg.data["ioda"]["tar_dir_template"]  # see note below
    tar_dir = template.format(
        base=exp.base_path,
        expid=exp.id,
        year=f"{dt.year:04d}",
        month=f"{dt.month:02d}",
    )
    if not os.path.isdir(tar_dir):
        raise FileNotFoundError(f"Resolved tar dir does not exist: {tar_dir!r}")
    return tar_dir




# ============================================================
# Experiment loading (with optional parallelism)
# ============================================================

def _load_all_experiments(cfg: Config, target: Target, start: datetime, end: datetime) -> Dict[str, TimeSeriesData]:
    """
    Build a TimeSeriesData per experiment for this target.
    """
    instrument = target.instrument
    varname = target.varname
    kx = target.kx
        
    ts_by_exp = {}
    for exp in cfg.experiments:
    #     log.info("Loading experiment %s (%s)", exp.id, exp.label)
        base_path = exp.base_path
        expid = exp.id
        file_type = exp.file_type
        if file_type == "ods":
            print("Loading ODS files...")
            dir_template = cfg.data.ods.dir_template
            file_pattern = cfg.data.ods.file_pattern
            file_time_format = cfg.data.ods.file_time_format
            ods_files = discover_ods_files(base_path,expid,instrument,dir_template,
                                           file_pattern,file_time_format,start,end)
            ts = build_ts_from_ods(ods_files,varname,kx,start,end)
            ts_by_exp[exp.id] = ts
            print("ODS files loaded")
        elif file_type == "ioda":
            print("Loading IODA files...")
            dir_template = cfg.data.ioda.dir_template
            file_pattern = cfg.data.ioda.file_pattern
            file_time_format = cfg.data.ioda.file_time_format
            ioda_files = discover_ioda_files(base_path,expid,instrument,dir_template,
                                           file_pattern,file_time_format,start,end)
            ts = build_ts_from_tar(ioda_files,instrument,varname,kx,start,end)
            ts.pass_data[0].data
            ts_by_exp[exp.id] = ts
            print("IODA files loaded")
        
    #     log.info("  loaded %d synoptic times for %s", len(ts.datetimes), exp.id)
    return ts_by_exp



def _aggregate_experiments(cfg, ts_by_exp: Dict[str, TimeSeriesData], 
                           start: datetime, end: datetime):
    ag_by_exp = {}
    for exp in ts_by_exp:
        ag_list = [aggregate_stats(ts_by_exp[exp], start, end)]
        pass_bin = aggregate_pass_binned(ts_by_exp[exp], start, end)
        pass_bin.data = replace(pass_bin.data, exp = exp)
        fail_bin = aggregate_fail_binned(ts_by_exp[exp], start, end)
        

        ag_list.extend([pass_bin, fail_bin])
        ag_by_exp[exp] = ag_list
    return ag_by_exp
    ...

# ============================================================
# Plot dispatch
# ============================================================

def _make_enabled_plots(cfg: Config, target: Target,
                        ts_by_exp: Dict[str, TimeSeriesData],ag_by_exp,
                        start: datetime, end: datetime) -> None:
    plots = cfg.plots

    

    if plots.statistics.enabled:
        _dispatch_statistics(cfg, ag_by_exp)


    if plots.coverage_map.enabled:
        _dispatch_coverage(cfg,ts_by_exp)

    if plots.comparison.enabled:
        _dispatch_comparison(cfg, ag_by_exp)

    if plots.time_series.enabled:
        _dispatch_time_series(cfg,target,ts_by_exp, start, end)


        ...





def _dispatch_statistics(cfg, ag_by_exp) -> None:
    iter_values = iter(ag_by_exp.values())
    ctl_values = next(iter_values)
    stats = ctl_values[0]
    pass_data = ctl_values[1]
    fail_data = ctl_values[2]
    fig = plot_stats(pass_data, fail_data, stats)
    ...


def _dispatch_time_series(cfg, target: Target, ts_by_exp, start, end) -> None:
    level = target.level
    for exp in ts_by_exp:
        ts = ts_by_exp[exp]
        fig = plot_series(ts, channel=level, start=start, end=end)
    ...

def _dispatch_coverage(cfg: Config,target:Target, ts_by_exp) -> None:
    level = target.level
    for exp in ts_by_exp:
        pass_data = ts_by_exp[exp].pass_data[0]
        fail_data = ts_by_exp[exp].fail_data[0]
    fig = plot_coverage(pass_data, fail_data, map_channel=level)

def _dispatch_comparison(cfg: Config, ag_by_exp) -> None:
    stats_type = cfg.plots.comparison.stat_type
    iter_values = iter(ag_by_exp.values())
    ctl_values = next(iter_values)
    ctl_stats = ctl_values[0]
    ctl_bins = ctl_values[1]
    exp_stats = next(iter_values)[0]

    fig = plot_purple(ctl_bins,ctl_stats,exp_stats, stats_type, cfg)


# ============================================================
# Output handling
# ============================================================

def _finalize(cfg: Config, fig, plot: str, target: Target, exp: str) -> None:
    """Honor output settings: save, show, or both."""
    out = cfg.output
    mode = out.get("mode", "show")

    if mode in ("save", "both"):
        os.makedirs(out["directory"], exist_ok=True)
        fname = out.get("filename_template", "{plot}_{instrument}_{exp}").format(
            plot=plot, instrument=target.instrument, exp=exp,
        )
        path = os.path.join(out["directory"], f"{fname}.{out.get('format','png')}")
        fig.savefig(path, dpi=out.get("dpi", 300), bbox_inches="tight")
        log.info("Saved %s", path)

    if mode in ("show", "both"):
        plt.show()

    plt.close(fig)   # free memory, important in long batch runs



