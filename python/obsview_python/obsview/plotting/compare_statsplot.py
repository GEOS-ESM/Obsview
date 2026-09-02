#Module containing functions to create 4 panel statistics plots comparing two different file sources
#This plot has the following capabilities:

    # 1. Compare two GSI or two JEDI output files
    #TODO 2. Compare single GSI output with single JEDI output
    #TODO 3. Compare HofX of GSI and JEDI from single IODA file
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

from ..processing.binning import BinnedData
from ..stats.statisticsdata import StatisticsData
from ..stats.calc_stats import count_obs_per_bin




def plot_compare_stats(pass_data1: BinnedData, stats1: StatisticsData, stats2: BinnedData):
    
    lev_type = pass_data1.data.lev_type

    fig = plt.figure(figsize=(10, 7))
  

    plt.subplot(2, 2, 1)
    _panel_nobs(pass_data1, stats1, stats2, lev_type)

    plt.subplot(2, 2, 2)
    _panel_resstats(pass_data1, stats1, stats2, lev_type)

    plt.subplot(2, 2, 3)
    _panel_jo(pass_data1, stats1, stats2, lev_type)

    plt.subplot(2, 2, 4)
    _panel_sigo(pass_data1, stats1, stats2, lev_type)
    #plt.subplots_adjust(top = 0.82)

    return fig


def _panel_nobs(pass_data1: BinnedData, stats1: StatisticsData, stats2: StatisticsData, lev_type):
    bin_centers = pass_data1.bin_centers
    bin_heights = pass_data1.bin_heights
    bar_width = bin_heights * 0.8   # single overlapping bar per level


    for i in range(len(bin_centers)):       #For each level...
        y = bin_centers[i]
        # Experiment 2 underneath.
        plt.barh(y, stats2.nobs[i], height=bar_width[i],
                 color="darkgreen", label="Experiment 2" if i == 0 else "", zorder=1)
        # Experiment 1 on top.
        plt.barh(y, stats1.nobs[i], height=bar_width[i],
                 color="lightgreen", alpha=0.6,
                 label="Experiment 1" if i == 0 else "", zorder=2)

    _apply_level_axis(plt.gca(), lev_type, "Observation Count")
    plt.tight_layout()

def _panel_resstats(binned_data, stats1, stats2, lev_type):
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels1 = ["Mean o-b: exp1", "Mean o-a: exp1", "RMS o-b: exp1", "RMS o-a: exp1"]
    labels2 = ["Mean o-b: exp2", "Mean o-a: exp2", "RMS o-b: exp2", "RMS o-a: exp2"]
    labeled1 = {"omb": False, "oma": False, "rms_omb": False, "rms_oma": False}
    labeled2 = {"omb": False, "oma": False, "rms_omb": False, "rms_oma": False}

    for i in range(len(bin_centers)):
        y = bin_centers[i]
        #Experiment 2
        plt.barh(y + offsets[0] * bar_width[i], stats2.mean_omb[i], height=bar_width[i],
                 color="cyan", label=labels2[0] if not labeled1["omb"] else "", zorder = 1)
        labeled1["omb"] = True
        plt.barh(y + offsets[2] * bar_width[i], stats2.mean_oma[i], height=bar_width[i],
                 color="orange", label=labels2[1] if not labeled1["oma"] else "", zorder = 1)
        labeled1["oma"] = True
        plt.barh(y + offsets[1] * bar_width[i], stats2.rms_omb[i], height=bar_width[i],
                 color="blue", label=labels2[2] if not labeled1["rms_omb"] else "", zorder = 1)
        labeled1["rms_omb"] = True
        plt.barh(y + offsets[3] * bar_width[i], stats2.rms_oma[i], height=bar_width[i],
                 color="red", label=labels2[3] if not labeled1["rms_oma"] else "", zorder = 1)
        labeled1["rms_oma"] = True
        
        
        #Experiment 1
        plt.barh(y + offsets[0] * bar_width[i], stats1.mean_omb[i], height=bar_width[i],
                 color="cyan",alpha=0.6, label=labels1[0] if not labeled2["omb"] else "", zorder = 2)
        labeled2["omb"] = True
        plt.barh(y + offsets[2] * bar_width[i], stats1.mean_oma[i], height=bar_width[i],
                 color="orange",alpha=0.6, label=labels1[1] if not labeled2["oma"] else "", zorder = 2)
        labeled2["oma"] = True
        plt.barh(y + offsets[1] * bar_width[i], stats1.rms_omb[i], height=bar_width[i],
                 color="blue",alpha=0.6, label=labels1[2] if not labeled2["rms_omb"] else "", zorder = 2)
        labeled2["rms_omb"] = True
        plt.barh(y + offsets[3] * bar_width[i], stats1.rms_oma[i], height=bar_width[i],
                 color="red",alpha=0.6, label=labels1[3] if not labeled2["rms_oma"] else "", zorder = 2)
        labeled2["rms_oma"] = True

    _apply_level_axis(plt.gca(), lev_type, "Mean & RMS of Obs Residuals")
    plt.tight_layout()


def _panel_jo(binned_data, stats1, stats2, lev_type):
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels1 = ["Jo(b)/p: exp1", "Jo(a)/p: exp1"]
    labels2 = ["Jo(b)/p: exp2", "Jo(a)/p: exp2"]
    labeled1 = {"job": False, "joa": False}
    labeled2 = {"job": False, "joa": False}

    for i in range(len(bin_centers)):
        y = bin_centers[i]
        
        #Experiment 2
        plt.barh(y + offsets[0] * bar_width[i], stats2.mean_job[i], height=bar_width[i],
                 color="blue", label=labels2[0] if not labeled2["job"] else "", zorder = 1)
        labeled2["job"] = True
        plt.barh(y + offsets[1] * bar_width[i], stats2.mean_joa[i], height=bar_width[i],
                 color="red", label=labels2[1] if not labeled2["joa"] else "", zorder = 1)
        labeled2["joa"] = True
        
        #Experiment 1
        plt.barh(y + offsets[0] * bar_width[i], stats1.mean_job[i], height=bar_width[i],
                 color="blue",alpha=0.6, label=labels1[0] if not labeled1["job"] else "", zorder = 2)
        labeled1["job"] = True
        plt.barh(y + offsets[1] * bar_width[i], stats1.mean_joa[i], height=bar_width[i],
                 color="red",alpha=0.6, label=labels1[1] if not labeled1["joa"] else "", zorder = 2)
        labeled1["joa"] = True

    _apply_level_axis(plt.gca(), lev_type, "Jo/p")
    plt.tight_layout()


def _panel_sigo(binned_data, stats1, stats2, lev_type):
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels1 = ["sigO: exp1", "esigO: exp1", "esigB: exp1"]
    labels2 = ["sigO: exp2", "esigO: exp2", "esigB: exp2"]
    labeled1 = {"sigo": False, "esigo": False, "esigb": False}
    labeled2 = {"sigo": False, "esigo": False, "esigb": False}

    for i in range(len(bin_centers)):
        y = bin_centers[i]
        #Experiment 2
        plt.barh(y + offsets[0] * bar_width[i], stats2.mean_sigo[i], height=bar_width[i],
                 color="cyan", label=labels2[0] if not labeled2["sigo"] else "", zorder = 1)
        labeled2["sigo"] = True
        plt.barh(y + offsets[1] * bar_width[i], stats2.mean_esigo[i], height=bar_width[i],
                 color="orange", label=labels2[1] if not labeled2["esigo"] else "", zorder = 1)
        labeled2["esigo"] = True
        plt.barh(y + offsets[2] * bar_width[i], stats2.mean_esigb[i], height=bar_width[i],
                 color="black", label=labels2[2] if not labeled2["esigb"] else "", zorder = 1)
        labeled2["esigb"] = True
        
        #Experiment 1
        plt.barh(y + offsets[0] * bar_width[i], stats1.mean_sigo[i], height=bar_width[i],
                 color="cyan",alpha=0.6, label=labels1[0] if not labeled1["sigo"] else "", zorder = 2)
        labeled1["sigo"] = True
        plt.barh(y + offsets[1] * bar_width[i], stats1.mean_esigo[i], height=bar_width[i],
                 color="orange",alpha=0.6, label=labels1[1] if not labeled1["esigo"] else "", zorder = 2)
        labeled1["esigo"] = True
        plt.barh(y + offsets[2] * bar_width[i], stats1.mean_esigb[i], height=bar_width[i],
                 color="black",alpha=0.6, label=labels1[2] if not labeled1["esigb"] else "", zorder = 2)
        labeled1["esigb"] = True

    _apply_level_axis(plt.gca(), lev_type, "Prescribed & Estimated Errors")
    plt.tight_layout()
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 10))



def _apply_level_axis(ax, lev_type: str, title_base: str):
    """
    Apply y-axis scaling, direction, label, and title suffix based on the
    vertical level type.

    - "pressure": log y-scale, inverted (1000 hPa at bottom), 'Pressure (hPa)'.
    - "channel" : linear y-scale, 'Channel'.
    """
    if lev_type == "pressure":
        ax.set_yscale("log")
        ax.invert_yaxis()                 # 1000 hPa at bottom, 0.1 hPa at top
        ax.set_ylabel("Pressure (hPa)")
        ax.set_title(f"{title_base} vs Pressure")
    elif lev_type == "channel":
        ax.set_ylabel("Channel")
        ax.set_title(f"{title_base} vs Channel")
    else:
        raise ValueError(f"Unknown lev_type: {lev_type!r}")

    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()    