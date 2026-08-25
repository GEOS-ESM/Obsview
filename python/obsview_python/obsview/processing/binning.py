#Module containing BinnedData object as well as functions to create bins and sort data by level
import numpy as np
from typing import Optional, List
from dataclasses import dataclass
from ..loading.observationdata import ObservationData
from .filtering import apply_filter

@dataclass
class BinnedData:
    data: ObservationData           #Rearranged data by bin
    bin_centers: np.ndarray         #Averaged value between bin levels (used for plotting)
           
    bin_labels: np.ndarray          #Array of each bin level (unique)
    bin_heights: np.ndarray
    bin_indices: Optional[np.ndarray] = None   
    #This should eventually replace the large ObservationData attribute
    nobs: Optional[np.ndarray] = None       #array with number of observations (int) for each vertical level 
    is_ts: Optional[bool] = None         #For plotting, tells us if the binned data came from a time series object
    ts_range: Optional[List[object]] = None   #For plotting, gives start and end datetimes










def create_bins(data: ObservationData) -> BinnedData:
    if data.lev_type == "channel":
        return create_channel_bins(data)
    elif data.lev_type == "pressure":
        return create_pressure_bins(data)
    raise ValueError(f"Unknown lev_type: {data.lev_type!r}")



def create_pressure_bins(data: ObservationData) -> BinnedData:
    NUM_BINS = 18                 # Hardcoded for now
    LEVLIM = [1000.0, 0.1]        # [bottom, top] in hPa, hardcoded for now
    if data.file_type == 'ioda':
        pressure_hpa = data.lev / 100.0         #IODA has different units than ODS
    else: 
        pressure_hpa = data.lev

    bins = np.logspace(np.log10(LEVLIM[1]), np.log10(LEVLIM[0]), num=NUM_BINS)

    raw_indices = np.digitize(pressure_hpa, bins)

    #Sort observations by pressure
    sort_indices = np.argsort(pressure_hpa)
    #Bin data
    binned_data = apply_filter(data, sort_indices)
    indices = raw_indices[sort_indices]

    # Bin centers/heights, one per interior bin (len(bins)-1 bins).
    centers = (bins[:-1] + bins[1:]) / 2
    heights = np.diff(bins)
    labels = centers

    obj = BinnedData(
        data=binned_data,
        bin_centers=centers,
        bin_indices=indices,
        bin_labels=labels,
        bin_heights=heights,
    )
    return obj

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