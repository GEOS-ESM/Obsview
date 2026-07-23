#Module containing BinnedData object as well as functions to create bins and sort data by level
import numpy as np
from dataclasses import dataclass
from ..loading.observationdata import ObservationData
from .filtering import apply_filter

@dataclass
class BinnedData:
    data: ObservationData           #Rearranged data by bin
    bin_centers: np.ndarray         #Averaged value between bin levels (used for plotting)
    bin_indices: np.ndarray         
    bin_labels: np.ndarray          #Array of each bin level (unique)
    bin_heights: np.ndarray
    #level_type: str (pressure or channel)

#TODO: define this function    
def create_pressure_bins():
    ...

def create_channel_bins(data: ObservationData) -> BinnedData:
    channels = np.unique(data.all_lev)   # Same as bin_labels
    
    # Create a sorting index that arranges data by channel
    sort_indices = np.argsort(data.lev)
    
    # Apply the sorting to get data organized by bin
    binned_data = apply_filter(data, sort_indices)
    
    # Now calculate bin indices based on the sorted data
    indices = np.searchsorted(channels, binned_data.lev)
    
    n = len(channels)
    centers = np.arange(1, n + 1)            
    heights = 0.95 * np.ones(n)      

    obj = BinnedData(
        data=binned_data,
        bin_centers=centers,
        bin_indices=indices,
        bin_labels=channels,
        bin_heights=heights
    )
    return obj