#Module containing plotting functions to create the 4 pannel statistics plot
import numpy as np
import matplotlib.pyplot as plt

from ..processing.binning import BinnedData
from ..stats.statisticsdata import StatisticsData
from ..stats.calc_stats import count_obs_per_bin

#TODO: Change y labels to reflect lev_type

MAIN_TITLE = "AMSU-A NOAA 19  |  2023-08-01 00:00 UTC"
def plot_stats(pass_data: BinnedData, fail_data:BinnedData, stats: StatisticsData, title: str = MAIN_TITLE):
    
    fig = plt.figure(figsize=(10, 7))  # hardcoded size for now
    fig.suptitle(title, fontsize = 14, fontweight = "bold")
    
    plt.subplot(2, 2, 1)
    _panel_nobs(pass_data, fail_data)

    plt.subplot(2, 2, 2)
    _panel_resstats(pass_data, stats)

    plt.subplot(2, 2, 3)
    _panel_jo(pass_data, stats)

    plt.subplot(2, 2, 4)
    _panel_sigo(pass_data, stats)

    return fig


def _panel_nobs(pass_data: BinnedData, fail_data: BinnedData):
    
    pass_nobs = count_obs_per_bin(pass_data)
    fail_nobs = count_obs_per_bin(fail_data)

    bin_centers = pass_data.bin_centers
    bin_heights = pass_data.bin_heights
    bar_width = bin_heights * 0.8   # wider single bar since we overlap now

    for i in range(len(bin_centers)):
        y = bin_centers[i]

        # Draw the "not used" (fail) bar first, fully opaque.
        plt.barh(
            y, fail_nobs[i],
            height=bar_width[i],
            color="red",
            label="not used" if i == 0 else "",
            zorder=1,
        )

        # Draw the "used" (pass) bar on top, at the SAME y, with opacity
        # so the red underneath is still visible.
        plt.barh(
            y, pass_nobs[i],
            height=bar_width[i],
            color="green",
            alpha=0.6,
            label="used" if i == 0 else "",
            zorder=2,
        )

    plt.ylabel("Channel")
    plt.title("Observation Count vs Channel")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()


def _panel_resstats(binned_data: BinnedData, stats: StatisticsData):
    """Panel 2: Mean & RMS of Obs Residuals vs Channel."""
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels = ["Mean o-b", "Mean o-a", "RMS o-b", "RMS o-a"]
    labeled = {"omb": False, "oma": False, "rms_omb": False, "rms_oma": False}

    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], stats.mean_omb[i], height=bar_width[i],
                 color="cyan", label=labels[0] if not labeled["omb"] else "")
        labeled["omb"] = True

        plt.barh(y + offsets[2] * bar_width[i], stats.mean_oma[i], height=bar_width[i],
                 color="orange", label=labels[1] if not labeled["oma"] else "")
        labeled["oma"] = True

        plt.barh(y + offsets[1] * bar_width[i], stats.rms_omb[i], height=bar_width[i],
                 color="blue", label=labels[2] if not labeled["rms_omb"] else "")
        labeled["rms_omb"] = True

        plt.barh(y + offsets[3] * bar_width[i], stats.rms_oma[i], height=bar_width[i],
                 color="red", label=labels[3] if not labeled["rms_oma"] else "")
        labeled["rms_oma"] = True

    plt.ylabel("Channel")
    plt.title("Mean & RMS of Obs Residuals vs Channel")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()


def _panel_jo(binned_data: BinnedData, stats: StatisticsData):
    """Panel 3: Jo/p vs Channel."""
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels = ["Jo(b)/p", "Jo(a)/p"]
    labeled = {"job": False, "joa": False}

    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], stats.mean_job[i], height=bar_width[i],
                 color="blue", label=labels[0] if not labeled["job"] else "")
        labeled["job"] = True

        plt.barh(y + offsets[1] * bar_width[i], stats.mean_joa[i], height=bar_width[i],
                 color="red", label=labels[1] if not labeled["joa"] else "")
        labeled["joa"] = True

    plt.ylabel("Channel")
    plt.title("Jo/p vs Channel")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()


def _panel_sigo(binned_data: BinnedData, stats: StatisticsData):
    """Panel 4: Prescribed & Estimated Errors vs Channel."""
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels = ["sigO", "esigO", "esigB"]
    labeled = {"sigo": False, "esigo": False, "esigb": False}

    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], stats.mean_sigo[i], height=bar_width[i],
                 color="cyan", label=labels[0] if not labeled["sigo"] else "")
        labeled["sigo"] = True

        plt.barh(y + offsets[1] * bar_width[i], stats.mean_esigo[i], height=bar_width[i],
                 color="orange", label=labels[1] if not labeled["esigo"] else "")
        labeled["esigo"] = True

        plt.barh(y + offsets[2] * bar_width[i], stats.mean_esigb[i], height=bar_width[i],
                 color="black", label=labels[2] if not labeled["esigb"] else "")
        labeled["esigb"] = True

    plt.ylabel("Channel")
    plt.title("Prescribed & Estimated Errors vs Channel")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 10))
    plt.legend()
