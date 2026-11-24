from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
from ..reference import Reference, Observation
from ..processing.binning import BinnedData
from ..stats.statisticsdata import StatisticsData
from ..config import Config

@dataclass(frozen=True)
class StatStyle:
    label: str
    color: str

# 1. Lookup tables for styling and human-readable names
STAT_STYLES = {
    "rms": StatStyle(label="RMS", color="blue"),
    "std": StatStyle(label="STD", color="purple"),
}

VAR_LABELS = {
    "omb": "O-B",
    "oma": "O-A",
    "obs": "Obs",
}

def plot_purple(
    ctl_bins: BinnedData,
    ctl_stats: StatisticsData,
    exp_stats: StatisticsData,
    stat_type: str,
    plot_var: str,
    cfg: Config,
    obs: Observation,
    ref: Reference,
):
    stat_key = stat_type.lower()
    var_key = plot_var.lower()

    # 2. Validate inputs against lookup dictionaries
    if stat_key not in STAT_STYLES:
        raise ValueError(f"Invalid stat_type '{stat_type}'. Expected one of {list(STAT_STYLES)}.")
    if var_key not in VAR_LABELS:
        raise ValueError(f"Invalid plot_var '{plot_var}'. Expected one of {list(VAR_LABELS)}.")

    style = STAT_STYLES[stat_key]
    var_label = VAR_LABELS[var_key]

    # 3. Dynamically fetch data attributes from StatisticsData
    ctl_mean = getattr(ctl_stats, f"mean_{var_key}")
    exp_mean = getattr(exp_stats, f"mean_{var_key}")
    diff = np.abs(exp_mean) - np.abs(ctl_mean)

    stat_attr = f"{stat_key}_{var_key}"
    ctl_spread = getattr(ctl_stats, stat_attr)
    exp_spread = getattr(exp_stats, stat_attr)
    spread_ratio = exp_spread / ctl_spread

    # --- Plotting Setup ---
    fig, ax1 = plt.subplots(figsize=(7, 8))
    y_pos = np.arange(1, len(ctl_bins.bin_labels) + 1)

    # Bottom X-axis (Mean Difference)
    diff_bound = cfg.plots.comparison.diff_bounds
    ax1.set_xlim(-diff_bound, diff_bound)
    ax1.barh(y_pos, diff, height=0.6, color="green", edgecolor="none", align="center", alpha=0.5)
    ax1.set_ylim(0, len(ctl_bins.bin_labels) + 1)
    ax1.axvline(0, color="purple", linestyle="-", linewidth=1, alpha=0.7)
    ax1.set_ylabel("Channel number", fontsize=12)
    ax1.set_xlabel(f"Mean {var_label} difference (E - C)\n<- Improvement | Deterioration ->", fontsize=11)

    labels = np.int16(ctl_bins.bin_labels)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(labels)

    if obs.level_type == "pressure":
        ax1.invert_yaxis()
        ax1.set_ylabel("Pressure level (hPa)", fontsize=12)

    ax1.tick_params(axis="x", which="major", labelsize=11, labelcolor="green")
    ax1.spines["top"].set_visible(True)
    ax1.spines["right"].set_visible(True)

    # Top X-axis (Ratio)
    ax2 = ax1.twiny()
    span = cfg.plots.comparison.ratio_bounds
    ax2.set_xlim(1 - span, span + 1)
    ax2.barh(
        y_pos,
        spread_ratio - 1,
        left=1,
        height=0.4,
        color=style.color,
        align="center",
        edgecolor="none",
        alpha=0.75,
    )
    ax2.set_xlabel(f"{style.label} {var_label} Ratio (E/C)", color=style.color, fontsize=11)
    ax2.tick_params(axis="x", labelcolor=style.color)

    title2 = f"Comparison of mean and {stat_key}: {var_label}"
    fig.text(0.5, 0.87, title2, fontsize=12, fontweight="bold", color=style.color, ha="center")

    # Horizontal grid marks
    for pos in y_pos:
        ax1.axhline(pos, color="gray", linestyle=":", linewidth=0.5, alpha=0.5)

    # Title
    data_src = ref.get_instrument_title(obs.instrument)
    data_type = ref.get_variable_title(obs.varname)
    ctl_name = cfg.experiments[0].id
    exp_name = cfg.experiments[1].id
    starttime, endtime = ctl_bins.ts_range[0], ctl_bins.ts_range[1]

    title = (
        f"{ctl_name} (ctl) vs. {exp_name} (exp)\n"
        f"{starttime.strftime('%d%b%Y %HZ')} - {endtime.strftime('%d%b%Y %HZ')}\n"
        f"{data_src} - {data_type}"
    )
    plt.suptitle(title, fontsize=14, fontweight="bold", y=0.98, ha="center")
    plt.tight_layout()
    plt.subplots_adjust(top=0.78)

    return fig