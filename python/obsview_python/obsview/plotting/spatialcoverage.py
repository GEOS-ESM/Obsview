#This module contains functions to create a spatial coverage (map) plot
import yaml
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path

from ..processing.binning import BinnedData


# ---- Easily-changeable coverage settings (hardcoded for now) ----
# Actual channel number to show on the MAP panel. Set to None to plot all
# channels. (The count/bar panel always shows every channel regardless.)
COVERAGE_MAP_LEVEL = 7

rc_path = Path(__file__).parent / '../resource.yaml'
with open(rc_path, "r") as f:
        rc = yaml.safe_load(f) 


def plot_coverage(pass_data: BinnedData,
                  fail_data: BinnedData,
                  map_level: int = COVERAGE_MAP_LEVEL):
   
    data_src = rc['KX_TO_DATASRC'].get(pass_data.data.kx)
    data_type =  rc['KT_TO_DATATYPE'].get(pass_data.data.kt)
    region = "Global"
    fig = plt.figure(figsize=(12, 10))

    time = pass_data.data.datetime
    title = f"Exp: {pass_data.data.exp} | {time.strftime('%d%b%Y %HZ')}"
    if pass_data.data.lev_type == "channel":
        title2 = f"{data_src} - {data_type}: channel {map_level} ({region})"
    else:
        title2 = f"{data_src} - {data_type}: level {map_level} hPa ({region})"


    plt.suptitle(title, fontsize=16, fontweight="bold",y = 0.97, ha="center")
    fig.text(0.5,0.91, title2, fontsize = 14, fontweight = "bold", color = "blue", ha = "center")
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1])

    # Map panel needs a cartopy projection axis.
    ax_map = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())
    _coverage_map_panel(ax_map, pass_data, fail_data, map_level)

    # Bar panel is a normal axis; always shows all levels.
    ax_bar = fig.add_subplot(gs[1])
    _coverage_count_panel(ax_bar, pass_data, fail_data)

    # Leave headroom for the suptitle so it doesn't overlap the map.
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def _coverage_map_panel(ax, pass_data: BinnedData, fail_data: BinnedData,
                        map_level: int):
    """
    Top panel: scatter observation locations on a world map.
    QC-fail (unused) drawn first in red; QC-pass (used) on top in green.
    If map_channel is not None, only observations on that channel are shown.
    """
    # Select observations for the requested channel (or keep all).
    if map_level is None:
        pass_lon, pass_lat = pass_data.data.lon, pass_data.data.lat
        fail_lon, fail_lat = fail_data.data.lon, fail_data.data.lat
        lev_str = "all channels"
    else:
        level_index = np.where(pass_data.bin_labels == map_level)[0]
        pass_sel = pass_data.bin_indices == level_index
        fail_sel = fail_data.bin_indices == level_index
        
        # pass_sel = pass_data.data.lev == map_level
        # fail_sel = fail_data.data.lev == map_level
        pass_lon, pass_lat = pass_data.data.lon[pass_sel], pass_data.data.lat[pass_sel]
        fail_lon, fail_lat = fail_data.data.lon[fail_sel], fail_data.data.lat[fail_sel]
        if pass_data.data.lev_type == "channel":
            lev_str = f"Channel #{int(map_level)}"
        else:
            lev_str = f"Level {int(map_level)} hPa"

    ax.set_global()

    # Base map features (mirrors radmap.py styling).
    ax.set_title(f"Observation Locations ({lev_str}) - QC-pass vs. QC-fail")
    ax.coastlines()
    ax.add_feature(cfeature.LAND, edgecolor="black", facecolor="lightgray")
    ax.add_feature(cfeature.OCEAN, facecolor="lightblue")
    ax.gridlines(draw_labels=True, linestyle="--", alpha=0.5)

    # Unused (QC-fail) first, semi-transparent, underneath.
    ax.scatter(
        fail_lon, fail_lat,
        s=5, color="red", alpha=0.3, label="not used",
        transform=ccrs.PlateCarree(), zorder=2,
    )

    # Used (QC-pass) on top, more opaque.
    ax.scatter(
        pass_lon, pass_lat,
        s=5, color="green", alpha=0.7, label="used",
        transform=ccrs.PlateCarree(), zorder=0,
    )

    ax.legend(loc="upper right")


def _coverage_count_panel(ax, pass_data: BinnedData, fail_data: BinnedData):
    """
    Bottom panel: per-channel count of used (green) vs. unused (red) obs,
    drawn as side-by-side vertical bars (like radmap.py's bar chart).
    Always shows every channel.
    """
    n_bins = len(pass_data.bin_labels)
    used_counts = np.bincount(pass_data.bin_indices, minlength=n_bins)
    unused_counts = np.bincount(fail_data.bin_indices, minlength=n_bins)

    x = np.arange(n_bins)
    width = 0.35

    ax.bar(x - width / 2, used_counts, width, label="used", color="green")
    ax.bar(x + width / 2, unused_counts, width, label="not used", color="red")

    # Label ticks with true level numbers (bin_labels)
    x_ax = np.arange(len(pass_data.bin_labels))
    labels = np.int32(pass_data.bin_labels)
    ax.set_xticks(x_ax)
    ax.set_xticklabels(labels)

    lev_type = pass_data.data.lev_type
    if lev_type == "channel":
         level = "channel"
    elif lev_type == "pressure":
         level = "level"
    ax.set_ylabel("Count")
    ax.set_title(f"Count of used and not-used observations per {level}")
    ax.grid(True, axis="y")
    ax.legend()