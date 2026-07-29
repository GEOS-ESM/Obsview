#Module for creating the time series plot (usual time length is one month) for data averaged globally 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import calendar
from datetime import datetime, timezone, timedelta

from ..loading.timeseriesdata import TimeSeriesData

SERIES_CHANNEL = 12
SERIES_TITLE = "ATMS N20 brightness temperatures: channel 12 (Global)"

def _default_month_range(datetimes):
    """
    Compute the default plotting range: the full calendar month containing
    the earliest datetime, from day 1 00Z through the last 18Z synoptic slot.
    """
    first = min(datetimes)
    start = datetime(first.year, first.month, 1, 0, 0, tzinfo=first.tzinfo)

    # Last day of that month, final 6-hourly synoptic slot (18Z).
    last_day = calendar.monthrange(first.year, first.month)[1]
    end = datetime(first.year, first.month, last_day, 18, 0, tzinfo=first.tzinfo)
    return start, end




def _channel_index(ts: TimeSeriesData, channel: int) -> int:
    """
    Find the array index of `channel` within the shared channel axis.
    Channels are assumed identical across files, so we read them from the
    first file's binned pass data.
    """
    labels = ts.pass_data[0].bin_labels
    idx = np.where(labels == channel)[0]
    if idx.size == 0:
        raise ValueError(f"Channel {channel} not found in {labels}")
    return int(idx[0])

def plot_series(ts: TimeSeriesData,
                channel: int = SERIES_CHANNEL,
                title: str = SERIES_TITLE,
                start: datetime = None,
                end: datetime = None):
    """
    Build the 3-panel time series figure for a single channel.

    The x-axis spans the full requested time range (default: the calendar
    month containing the first file), so missing synoptic times appear as
    empty gaps rather than compressing the axis.

    Parameters
    ----------
    start, end : datetime or None
        Explicit plotting range. If None, defaults to the full calendar
        month containing the earliest datetime in `ts`.
    """
    ci = _channel_index(ts, channel)

    # Full time range for the x-axis (independent of how many files exist).
    if start is None or end is None:
        default_start, default_end = _default_month_range(ts.datetimes)
        start = start or default_start
        end = end or default_end

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    fig.suptitle(title, fontsize=13, fontweight="bold", color="blue")

    times = ts.datetimes

    _series_counts_panel(ax1, ts, times, ci)
    _series_residuals_panel(ax2, ts, times, ci)
    _series_cost_panel(ax3, ts, times, ci)

    # Pin every panel to the full range so gaps are preserved.
    # (sharex=True means setting one propagates to all three.)
    ax3.set_xlim(start, end)

    # Sparse, readable date ticks across the full month, matching the figure
    # style (odd days: 01Mar, 03Mar, ...).
    ax3.xaxis.set_major_locator(mdates.DayLocator(bymonthday=range(1, 32, 2)))
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%d%b"))
    fig.autofmt_xdate(rotation=0, ha="center")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


def _series_counts_panel(ax, ts: TimeSeriesData, times, ci: int):
    """Panel 1: used vs. not-used observation counts for the channel."""
    # Used counts come from the pass StatisticsData per time.
    used = np.array([s.nobs[ci] for s in ts.pass_stats])

    # Fail data has no StatisticsData, so count per time from bin_indices.
    unused = np.array([
        np.bincount(bd.bin_indices,
                    minlength=len(bd.bin_labels))[ci]
        for bd in ts.fail_data
    ])

    # 6-hour synoptic spacing = 0.25 day; narrow slightly so bars separate.
    w = 0.20

    # Not-used underneath (red), used on top (green).
    ax.bar(times, unused, width=w, color="red", label="Not used", zorder=1)
    ax.bar(times, used, width=w, color="green", alpha=0.9,
           label="Used", zorder=2)

    ax.set_title("Data counts:", loc="left", fontsize=10)
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    ax.legend(loc="upper right", fontsize=8, ncol=2)


def _series_residuals_panel(ax, ts: TimeSeriesData, times, ci: int):
    """Panel 2: rms & mean of O-B and O-A for the channel."""
    rms_omb  = np.array([s.rms_omb[ci]  for s in ts.pass_stats])
    rms_oma  = np.array([s.rms_oma[ci]  for s in ts.pass_stats])
    mean_omb = np.array([s.mean_omb[ci] for s in ts.pass_stats])
    mean_oma = np.array([s.mean_oma[ci] for s in ts.pass_stats])

    w = 0.20
    ax.bar(times, rms_omb,  width=w, color="blue",   label="rms(O-B)",  zorder=1)
    ax.bar(times, rms_oma,  width=w, color="red",    label="rms(O-A)",  zorder=2)
    ax.bar(times, mean_omb, width=w, color="cyan",   label="mean(O-B)", zorder=3)
    ax.bar(times, mean_oma, width=w, color="orange", label="mean(O-A)", zorder=4)

    ax.axhline(0, color="black", linewidth=0.6)
    ax.set_title("Data residuals:", loc="left", fontsize=10)
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    ax.legend(loc="upper right", fontsize=8, ncol=4)


def _series_cost_panel(ax, ts: TimeSeriesData, times, ci: int):
    """Panel 3: Jo/p for O-B and O-A for the channel."""
    job = np.array([s.mean_job[ci] for s in ts.pass_stats])
    joa = np.array([s.mean_joa[ci] for s in ts.pass_stats])

    w = 0.20
    ax.bar(times, job, width=w, color="blue", label="Jo(O-B)/p", zorder=1)
    ax.bar(times, joa, width=w, color="red",  label="Jo(O-A)/p", zorder=2)

    ax.set_title("Normalized cost:", loc="left", fontsize=10)
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)
    ax.legend(loc="upper right", fontsize=8, ncol=2)


