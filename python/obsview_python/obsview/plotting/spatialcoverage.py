#This module contains functions to create a spatial coverage (map) plot
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from ..processing.binning import BinnedData
from ..config import KX_TO_DATASRC, KT_TO_DATATYPE

# ---- Easily-changeable coverage settings (hardcoded for now) ----
# Actual channel number to show on the MAP panel. Set to None to plot all
# channels. (The count/bar panel always shows every channel regardless.)
COVERAGE_MAP_CHANNEL = 7




def plot_coverage(pass_data: BinnedData,
                  fail_data: BinnedData,
                  map_channel: int = COVERAGE_MAP_CHANNEL):
    """
    Build the 2-panel coverage figure from the QC-pass and QC-fail
    BinnedData objects.

    Panel 1 (top): map of observation locations (used vs. unused) for a
                   single channel, or all channels if map_channel is None.
    Panel 2 (bottom): per-channel count of used vs. unused observations
                      (always shows all channels).

    Parameters
    ----------
    map_channel : int or None
        Actual channel number to display on the map. None => all channels.
    title : str
        Figure super-title (instrument/satellite + date/time).

    Returns the matplotlib Figure so the caller can plt.show() or savefig().
    """
    data_src = KX_TO_DATASRC.get(pass_data.data.kx)
    data_type =  KT_TO_DATATYPE.get(pass_data.data.kt)
    region = "Global"
    fig = plt.figure(figsize=(12, 10))

    time = pass_data.data.datetime
    title = f"Exp: {pass_data.data.exp} | {time.strftime('%d%b%Y %HZ')}"
    title2 = f"{data_src} - {data_type}: channel: {map_channel} ({region})"


    plt.suptitle(title, fontsize=16, fontweight="bold",y = 0.97, ha="center")
    fig.text(0.5,0.91, title2, fontsize = 14, fontweight = "bold", color = "blue", ha = "center")
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1])

    # Map panel needs a cartopy projection axis.
    ax_map = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())
    _coverage_map_panel(ax_map, pass_data, fail_data, map_channel)

    # Bar panel is a normal axis; always shows all channels.
    ax_bar = fig.add_subplot(gs[1])
    _coverage_count_panel(ax_bar, pass_data, fail_data)

    # Leave headroom for the suptitle so it doesn't overlap the map.
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def _coverage_map_panel(ax, pass_data: BinnedData, fail_data: BinnedData,
                        map_channel: int):
    """
    Top panel: scatter observation locations on a world map.
    QC-fail (unused) drawn first in red; QC-pass (used) on top in green.
    If map_channel is not None, only observations on that channel are shown.
    """
    # Select observations for the requested channel (or keep all).
    if map_channel is None:
        pass_lon, pass_lat = pass_data.data.lon, pass_data.data.lat
        fail_lon, fail_lat = fail_data.data.lon, fail_data.data.lat
        chan_str = "all channels"
    else:
        pass_sel = pass_data.data.lev == map_channel
        fail_sel = fail_data.data.lev == map_channel
        pass_lon, pass_lat = pass_data.data.lon[pass_sel], pass_data.data.lat[pass_sel]
        fail_lon, fail_lat = fail_data.data.lon[fail_sel], fail_data.data.lat[fail_sel]
        chan_str = f"Channel #{int(map_channel)}"

    # Base map features (mirrors radmap.py styling).
    ax.set_title(f"Observation Locations ({chan_str}) - QC-pass vs. QC-fail")
    ax.coastlines()
    ax.add_feature(cfeature.LAND, edgecolor="black", facecolor="lightgray")
    ax.add_feature(cfeature.OCEAN, facecolor="lightblue")
    ax.gridlines(draw_labels=True, linestyle="--", alpha=0.5)

    # Unused (QC-fail) first, semi-transparent, underneath.
    ax.scatter(
        fail_lon, fail_lat,
        s=5, color="red", alpha=0.3, label="not used",
        transform=ccrs.PlateCarree(), zorder=1,
    )

    # Used (QC-pass) on top, more opaque.
    ax.scatter(
        pass_lon, pass_lat,
        s=5, color="green", alpha=0.7, label="used",
        transform=ccrs.PlateCarree(), zorder=3,
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

    # Label ticks with true channel numbers (bin_labels), like radmap.py.
    ax.set_xticks(x)
    ax.set_xticklabels([f"#{int(ch)}" for ch in pass_data.bin_labels])

    ax.set_ylabel("Count")
    ax.set_title("Count of used and not-used observations per Channel")
    ax.grid(True, axis="y")
    ax.legend()