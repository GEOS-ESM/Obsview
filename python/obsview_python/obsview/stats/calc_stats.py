#Module containing calculation functions for statistics
import numpy as np
from .statisticsdata import StatisticsData
from ..processing.binning import BinnedData

def count_obs_per_bin(binned_data: BinnedData) -> np.ndarray:
    """Return per-bin observation counts only (no derived-variable access)."""
    n_bins = len(binned_data.bin_labels)
    return np.bincount(binned_data.bin_indices, minlength=n_bins)



def calculate_stats(binned_data: BinnedData) -> StatisticsData:
    n_bins = len(binned_data.bin_labels)
    
    # Count observations per bin
    nobs = np.bincount(binned_data.bin_indices, minlength=n_bins)
    
    # Find the start index of each bin in the sorted data
    bin_starts = np.searchsorted(binned_data.bin_indices, np.arange(n_bins))
    bin_ends = np.append(bin_starts[1:], len(binned_data.bin_indices))
    
    # Pre-allocate arrays
    mean_omb = np.full(n_bins, np.nan)
    mean_oma = np.full(n_bins, np.nan)
    rms_omb = np.full(n_bins, np.nan)
    rms_oma = np.full(n_bins, np.nan)
    mean_job = np.full(n_bins, np.nan)
    mean_joa = np.full(n_bins, np.nan)
    mean_sigo = np.full(n_bins, np.nan)
    mean_esigo = np.full(n_bins, np.nan)
    mean_esigb = np.full(n_bins, np.nan)
    
    # Calculate statistics for each bin
    for i in range(n_bins):
        if nobs[i] > 0:
            start, end = bin_starts[i], bin_ends[i]
            
            mean_omb[i] = np.mean(binned_data.data.omb[start:end])
            mean_oma[i] = np.mean(binned_data.data.oma[start:end])
            rms_omb[i] = np.sqrt(np.mean(binned_data.data.omb[start:end]**2))
            rms_oma[i] = np.sqrt(np.mean(binned_data.data.oma[start:end]**2))
            mean_job[i] = np.mean(binned_data.data.job[start:end])
            mean_joa[i] = np.mean(binned_data.data.joa[start:end])
            mean_sigo[i] = np.mean(binned_data.data.sigo[start:end])
            mean_esigo[i] = np.sqrt(np.abs(np.mean(binned_data.data.esigo[start:end])))
            mean_esigb[i] = np.sqrt(np.abs(np.mean(binned_data.data.esigb[start:end])))
    
    return StatisticsData(
        nobs=nobs,
        mean_omb=mean_omb,
        mean_oma=mean_oma,
        rms_omb=rms_omb,
        rms_oma=rms_oma,
        mean_job=mean_job,
        mean_joa=mean_joa,
        mean_sigo=mean_sigo,
        mean_esigo=mean_esigo,
        mean_esigb=mean_esigb
    )
    ...
