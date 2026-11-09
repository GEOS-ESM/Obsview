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
import logging
import numpy as np
from datetime import datetime
from typing import Dict
from dataclasses import replace

import matplotlib.pyplot as plt

from .config import Config, Experiment, Target
from .reference import Reference, Observation

# Time series construction + plotting (existing modules).
from .utils import _parse_time, _dt_to_int
from .loading.timeseriesdata import TimeSeriesData
from .loading.search import discover_ods_files, discover_ioda_files
from .loading.buildtimeseries import build_ts_from_ods, build_ts_from_tar
from .stats.aggregate import aggregate_stats, aggregate_pass_binned, aggregate_fail_binned
from .plotting.statsplot import plot_stats
from .plotting.timeseries import plot_series
from .plotting.spatialcoverage import plot_coverage
from .plotting.purple import plot_purple


log = logging.getLogger("obsview")


# ============================================================
# Entry point
# ============================================================

def run(cfg: Config, ref: Reference) -> None:
    #_configure_logging(cfg)

    #log.info("Starting run: file_type=%s, %d experiment(s), %d target(s)",
             #cfg.file_type, len(cfg.experiments), len(cfg.targets))

    start = _parse_time(cfg.selection.time_range.start)
    end = _parse_time(cfg.selection.time_range.end)

    # Outer loop: each instrument/obtype the user asked to analyze.
    for target in cfg.selection.targets:

        # Load each experiment's time series for this target.
        ts_by_exp: Dict[str, TimeSeriesData] = _load_all_experiments(cfg, ref, target, start, end)
        ag_by_exp = _aggregate_experiments(cfg,ts_by_exp, start, end)


        # Dispatch to enabled plots
        _make_enabled_plots(cfg, target, ts_by_exp, ag_by_exp, start, end, ref)
        plt.show()

    #log.info("Run complete.")
    print(f"Run complete")

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
# Experiment loading (with optional parallelism)
# ============================================================

def _load_all_experiments(cfg: Config, ref: Reference, target: Target, start: datetime, end: datetime) -> Dict[str, TimeSeriesData]:
    """
    Build a TimeSeriesData per experiment for this target.
    """
    #Map config target to reference observation
    obs = ref.get_observation(target.instrument, target.varname)

    instrument = obs.instrument
    varname = obs.varname
    kx = obs.kx
        
    ts_by_exp = {}

    for exp in cfg.experiments:
        base_path = exp.base_path
        expid = exp.id
        file_type = exp.file_type
        if file_type == "ods":
            print(f"Loading ODS files for {exp.id}: {instrument}...")
            dir_template = cfg.data.ods.dir_template
            file_pattern = cfg.data.ods.file_pattern
            file_time_format = cfg.data.ods.file_time_format
            ods_files = discover_ods_files(base_path,expid,instrument,dir_template,
                                           file_pattern,file_time_format,start,end)
            ts = build_ts_from_ods(ods_files,varname,kx,start,end,obs,ref)
            ts_by_exp[exp.id] = ts
            print("ODS files loaded")
        elif file_type == "ioda":
            print(f"Loading IODA files for {exp.id}: {instrument}...")
            dir_template = cfg.data.ioda.dir_template
            file_pattern = cfg.data.ioda.file_pattern
            file_time_format = cfg.data.ioda.file_time_format
            ioda_files = discover_ioda_files(base_path,expid,instrument,dir_template,
                                           file_pattern,file_time_format,start,end)
            ts = build_ts_from_tar(ioda_files,instrument,varname,kx,start,end,obs,ref)
            ts.pass_data[0].data
            ts_by_exp[exp.id] = ts
            print("IODA files loaded")

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
                        start: datetime, end: datetime, ref: Reference) -> None:
    plots = cfg.plots

    

    if plots.statistics.enabled:
        _dispatch_statistics(cfg, target, ag_by_exp, start, end, ref)


    if plots.coverage_map.enabled:
        _dispatch_coverage(cfg,target,ts_by_exp, start, end, ref)

    if plots.comparison.enabled:
        _dispatch_comparison(cfg, ag_by_exp, target, ref)

    if plots.time_series.enabled:
        _dispatch_time_series(cfg,target,ts_by_exp, start, end, ref)


        ...





def _dispatch_statistics(cfg: Config, target: Target, ag_by_exp, start: datetime, end: datetime, ref: Reference) -> None:
    #TODO: Loop through experiments instead of picking the first one
    iter_values = iter(ag_by_exp.values())
    ctl_values = next(iter_values)
    stats = ctl_values[0]
    pass_data = ctl_values[1]
    fail_data = ctl_values[2]

    exp = next(iter(ag_by_exp))
    start = _dt_to_int(start)
    end = _dt_to_int(end)

    obs = ref.get_observation(target.instrument, target.varname)

    fig = plot_stats(pass_data, fail_data, stats, obs, ref)
    _finalize(cfg,fig,"stats",target,exp, start, end)
    ...


def _dispatch_time_series(cfg, target, ts_by_exp, start: datetime, end: datetime, ref: Reference) -> None:
    level = target.level
    start = _dt_to_int(start)
    end = _dt_to_int(end)

    obs = ref.get_observation(target.instrument, target.varname)

    for exp in ts_by_exp:
        ts = ts_by_exp[exp]
        if level == "all":
            for lev in ts.pass_data[0].bin_labels:
                level_num = np.where(ts.pass_data[0].bin_labels == lev)[0][0]+1
                fig = plot_series(ts, obs, ref, channel=lev, start=None, end=None)
                _finalize_lev(cfg, fig, "timeseries", target, exp, start, end, level_num)
        else:
            level_num = np.where(ts.pass_data[0].bin_labels == level)[0][0]+1
            fig = plot_series(ts, obs, ref, channel=level, start=None, end=None)
            _finalize_lev(cfg, fig, "timeseries", target, exp, start, end, level_num)


    ...

def _dispatch_coverage(cfg: Config,target:Target, ts_by_exp, start, end, ref: Reference) -> None:
    level = target.level
    start = _dt_to_int(start)
    end = _dt_to_int(end)

    obs = ref.get_observation(target.instrument, target.varname)


    for exp in ts_by_exp:
        pass_data = ts_by_exp[exp].pass_data[0]
        fail_data = ts_by_exp[exp].fail_data[0]
        if level == "all":
            for lev in pass_data.bin_labels:
                level_num = np.where(pass_data.bin_labels == lev)[0][0]+1
                fig = plot_coverage(pass_data, fail_data, obs, ref, map_level=lev)
                _finalize_lev(cfg, fig, "coverage", target, exp, start, end, level_num)
        else: 
            level_num = np.where(pass_data.bin_labels == level)[0][0]+1
            fig = plot_coverage(pass_data, fail_data, obs, ref, map_level=level)
            _finalize_lev(cfg, fig, "coverage", target, exp, start, end, level_num)


def _dispatch_comparison(cfg: Config, ag_by_exp, target: Target, ref: Reference) -> None:
    stats_type = cfg.plots.comparison.stat_type
    iter_values = iter(ag_by_exp.values())
    ctl_values = next(iter_values)
    ctl_stats = ctl_values[0]
    ctl_bins = ctl_values[1]
    exp_stats = next(iter_values)[0]

    obs = ref.get_observation(target.instrument, target.varname)

    fig = plot_purple(ctl_bins,ctl_stats,exp_stats, stats_type, cfg, obs, ref)


# ============================================================
# Output handling
# ============================================================

def _finalize(cfg: Config, fig, plot: str, target: Target, exp: str, start: int, end: int) -> None:
    """Honor output settings: save, show, or both."""
    out = cfg.output
    mode = out.mode

    if mode in ("save", "both"):
        os.makedirs(out.save_directory, exist_ok=True)
        fname = out.filename_template.format(
            plot=plot, instrument=target.instrument, exp=exp,start = start, end = end
        )
        path = os.path.join(out.save_directory, f"{fname}.{out.format}")
        fig.savefig(path, dpi=out.dpi, bbox_inches="tight")
        log.info("Saved %s", path)

    if mode in ("show", "both"):
        plt.show()
    plt.close(fig)

#Finalize function for plots that are level specific(e.g. time series)
def _finalize_lev(cfg: Config, fig, plot: str, target: Target, exp: str, start: int, end: int, level) -> None:
    """Honor output settings: save, show, or both."""
    out = cfg.output
    mode = out.mode

    if mode in ("save", "both"):
        os.makedirs(out.save_directory, exist_ok=True)
        fname = out.filename_lev_template.format(
            plot=plot, instrument=target.instrument, exp=exp,start = start, end = end, level = level
        )
        path = os.path.join(out.save_directory, f"{fname}.{out.format}")
        fig.savefig(path, dpi=out.dpi, bbox_inches="tight")
        log.info("Saved %s", path)

    if mode in ("show", "both"):
        plt.show()
    plt.close(fig)   



