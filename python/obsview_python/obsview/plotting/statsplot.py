#Module containing plotting functions to create the 4 pannel statistics plot
import yaml
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from datetime import datetime
from pathlib import Path

from ..processing.binning import BinnedData
from ..stats.statisticsdata import StatisticsData
from ..stats.calc_stats import count_obs_per_bin



rc_path = Path(__file__).parent / '../resource.yaml'
with open(rc_path, "r") as f:
        rc = yaml.safe_load(f) 

#Adjust bar height (0-1 relative to level size)
BAR_H = 0.8

def plot_stats(pass_data: BinnedData, fail_data: BinnedData, stats: StatisticsData):
    """
    4-panel statistics figure. Y-axis scaling/labels adapt to the vertical
    level type ('channel' or 'pressure') carried by the data.
    """
    lev_type = pass_data.data.lev_type

    data_src = rc['KX_TO_DATASRC'].get(pass_data.data.kx)
    data_type =  rc['KT_TO_DATATYPE'].get(pass_data.data.kt)

    #For plotting single time or time range
    if pass_data.is_ts:
        starttime = pass_data.ts_range[0]
        endtime = pass_data.ts_range[1]
        title = f"Exp: {pass_data.data.exp} | {starttime.strftime('%d%b%Y %HZ')} - {endtime.strftime('%d%b%Y %HZ')}"
    else:
        time = pass_data.data.datetime
        title = f"Exp: {pass_data.data.exp} | {time.strftime('%d%b%Y %HZ')}"
    
    title2 = f"{data_src} - {data_type}"


    fig = plt.figure(figsize=(10, 7))
    plt.suptitle(title, fontsize=16, fontweight="bold",y = 0.98, ha="center")
    fig.text(0.5,0.91, title2, fontsize = 16, fontweight = "bold", color = "blue", ha = "center")
    

    plt.subplot(2, 2, 1)
    _panel_nobs(pass_data, fail_data, stats, lev_type)

    plt.subplot(2, 2, 2)
    _panel_resstats(pass_data, stats, lev_type)

    plt.subplot(2, 2, 3)
    _panel_jo(pass_data, stats, lev_type)

    plt.subplot(2, 2, 4)
    _panel_sigo(pass_data, stats, lev_type)
    plt.subplots_adjust(top = 0.82)

    return fig


def _panel_nobs(pass_data: BinnedData, fail_data: BinnedData, stats: StatisticsData, lev_type):
    bar_width = BAR_H

    # Fail counts per bin (fail_data has no StatisticsData).
    n_bins = len(fail_data.bin_labels)
    fail_nobs = fail_data.nobs

    for i in range(n_bins):       #For each level...
        y = range(n_bins)[i]
        # Red (not used) underneath.
        plt.barh(y, fail_nobs[i], height=bar_width,
                 color="red", label="not used" if i == 0 else "", zorder=1)
        # Green (used) on top, semi-transparent so red shows through.
        plt.barh(y, stats.nobs[i], height=bar_width,
                 color="green", alpha=0.6,
                 label="used" if i == 0 else "", zorder=2)

    _apply_level_axis(plt.gca(), pass_data, lev_type, "Observation Count")
    plt.tight_layout()

def _panel_resstats(binned_data, stats, lev_type):
    n_bins = len(binned_data.bin_labels)
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = BAR_H/2
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels = ["Mean o-b", "Mean o-a", "RMS o-b", "RMS o-a"]
    labeled = {"omb": False, "oma": False, "rms_omb": False, "rms_oma": False}

    for i in range(n_bins):
        y = range(n_bins)[i]
        plt.barh(y + offsets[0] * bar_width, stats.mean_omb[i], height=bar_width,
                 color="cyan", label=labels[0] if not labeled["omb"] else "")
        labeled["omb"] = True
        plt.barh(y + offsets[2] * bar_width, stats.mean_oma[i], height=bar_width,
                 color="orange", label=labels[1] if not labeled["oma"] else "")
        labeled["oma"] = True
        plt.barh(y + offsets[1] * bar_width, stats.rms_omb[i], height=bar_width,
                 color="blue", label=labels[2] if not labeled["rms_omb"] else "")
        labeled["rms_omb"] = True
        plt.barh(y + offsets[3] * bar_width, stats.rms_oma[i], height=bar_width,
                 color="red", label=labels[3] if not labeled["rms_oma"] else "")
        labeled["rms_oma"] = True

    _apply_level_axis(plt.gca(), binned_data, lev_type, "Mean & RMS of Obs Residuals")
    plt.tight_layout()


def _panel_jo(binned_data, stats, lev_type):
    n_bins = len(binned_data.bin_labels)
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = BAR_H/2
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels = ["Jo(b)/p", "Jo(a)/p"]
    labeled = {"job": False, "joa": False}

    for i in range(n_bins):
        y = range(n_bins)[i]
        plt.barh(y + offsets[0] * bar_width, stats.mean_job[i], height=bar_width,
                 color="blue", label=labels[0] if not labeled["job"] else "")
        labeled["job"] = True
        plt.barh(y + offsets[1] * bar_width, stats.mean_joa[i], height=bar_width,
                 color="red", label=labels[1] if not labeled["joa"] else "")
        labeled["joa"] = True

    _apply_level_axis(plt.gca(),binned_data, lev_type, "Jo/p")
    plt.tight_layout()


def _panel_sigo(binned_data, stats, lev_type):
    n_bins = len(binned_data.bin_labels)
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = BAR_H/2
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels = ["sigO", "esigO", "esigB"]
    labeled = {"sigo": False, "esigo": False, "esigb": False}

    for i in range(n_bins):
        y = range(n_bins)[i]
        plt.barh(y + offsets[0] * bar_width, stats.mean_sigo[i], height=bar_width,
                 color="cyan", label=labels[0] if not labeled["sigo"] else "")
        labeled["sigo"] = True
        plt.barh(y + offsets[1] * bar_width, stats.mean_esigo[i], height=bar_width,
                 color="orange", label=labels[1] if not labeled["esigo"] else "")
        labeled["esigo"] = True
        plt.barh(y + offsets[2] * bar_width, stats.mean_esigb[i], height=bar_width,
                 color="black", label=labels[2] if not labeled["esigb"] else "")
        labeled["esigb"] = True

    _apply_level_axis(plt.gca(),binned_data, lev_type, "Prescribed & Estimated Errors")
    plt.tight_layout()
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 10))



def _apply_level_axis(ax, pass_data: BinnedData, lev_type: str, title_base: str):
    y_ax = np.arange(len(pass_data.bin_labels))
    labels = np.int32(pass_data.bin_labels)
    ax.set_yticks(y_ax)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.5, len(y_ax) - 0.5)
    if lev_type == "pressure":
        ax.invert_yaxis()                 
        ax.set_ylabel("Pressure (hPa)")
        ax.set_title(f"{title_base} vs Pressure")
    elif lev_type == "channel":
        ax.set_ylabel("Channel")
        ax.set_title(f"{title_base} vs Channel")
    else:
        raise ValueError(f"Unknown lev_type: {lev_type!r}")

    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()