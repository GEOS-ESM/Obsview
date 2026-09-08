#Module containing plotting functions similar to the "purple" plots
#Plots time averaged statistics
import numpy as np
import matplotlib.pyplot as plt

from ..processing.binning import BinnedData
from ..stats.statisticsdata import StatisticsData






def plot_purple(ctl_bins: BinnedData, ctl_stats: StatisticsData, exp_stats:StatisticsData):
    fig, ax1 = plt.subplots(figsize=(7, 8))

    ctl = ctl_stats.mean_omb
    exp = exp_stats.mean_omb
    diff = (np.abs(exp) - np.abs(ctl))
    channels = ctl_bins.bin_labels


    rms_ratio = (exp_stats.rms_omb/ctl_stats.rms_omb)
    # Mock Confidence Intervals
    lower_ci = rms_ratio - 0.02
    upper_ci = rms_ratio + 0.03

    # Calculate error lengths (distances from the center dot)
    lower_error = rms_ratio - lower_ci
    upper_error = upper_ci - rms_ratio
    error_spans = np.array([lower_error, upper_error])


    bar_colors = ['green' if val < 0 else 'red' for val in diff]
    ratio_bar_colors = ['green' if val < 1 else 'red' for val in rms_ratio]
    max_val = np.nanmax(np.abs(diff))

# 2. Add an optional padding margin (e.g., 10% extra space) so bars don't touch the edges
    padding = 1.10 
    x_limit = max_val * padding

# 3. Apply symmetric limits to keep zero dead-center
    ax1.set_xlim(-x_limit, x_limit)

    #Plot difference ratio
    ax1.barh(channels, diff, height=0.6, color="green", edgecolor='none', align='center', alpha = 0.5)
    
    ax1.set_ylim(0, np.max(channels)+1)
    ax1.axvline(0, color='purple', linestyle='-', linewidth=1, alpha=0.7) # Reference line at zero


    ax1.set_ylabel('Channel number', fontsize=12)
    ax1.set_xlabel(f'Mean O-B difference (E - C)\n<- Improvement | Deterioration ->', fontsize=11)

    
    ax1.set_yticks(np.arange(2, 24, 2))
    

    #Styling
    ax1.tick_params(axis='x', which='major', labelsize=11, labelcolor = "green")
    ax1.spines['top'].set_visible(True)
    ax1.spines['right'].set_visible(True)


    # 3. Create Twin Axis for RMS (Top X-axis)
    ax2 = ax1.twiny()  

    #ax2.errorbar(rms_ratio, channels, xerr=error_spans, fmt='none', ecolor='blue',elinewidth=1.2,capsize=5, alpha=0.6)
    #ax2.plot(rms_ratio, channels, color='blue', linestyle='', marker='o', linewidth=1.5, label='RMS Diff')

    ax2.barh(channels, rms_ratio-1,left = 1, height=0.4, color="blue",align = "center", edgecolor='none', alpha = 0.75)
    
    span = 0.15
    ax2.set_xlim(1-span, span+1)  # Scales from 0 to slightly above max RMS

# Labels for Axis 2 (Top)
    ax2.set_xlabel('RMS O-B Ratio (E/C)', color='blue', fontsize=11)
    ax2.tick_params(axis='x', labelcolor='blue')


    #Horizontal marks:
    for channel in channels:
        ax1.axhline(channel, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)

    #Title
    title = f"x0053RPY (ctl) vs. x0054 (exp)\nComparison of mean and rms: O-B"
    plt.suptitle(title, fontsize = 14, fontweight = "bold", ha = "center")
    plt.tight_layout()
    ...



