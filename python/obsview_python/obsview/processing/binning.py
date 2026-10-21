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

DEFAULT_PRESSURE_BINS = np.array([5.0,7.0,10.0,20.0,30.0,50.0,70.0,100.0,150.0,200.0,250.0,300.0,400.0,500.0,700.0,850.0,925.0,1000.0])



def create_pressure_bins(data: ObservationData, bins = DEFAULT_PRESSURE_BINS) -> BinnedData:
    if data.file_type == 'ioda':
        pressure_hpa = data.lev / 100.0
    else: 
        pressure_hpa = data.lev

    # 1. Calculate midpoints between standard levels to create edges
    midpoints = (bins[:-1] + bins[1:]) / 2
    
    # 2. Add outer padding: 0.0 hPa at the top of atmosphere, np.inf at the surface
    # This creates exactly n+1 edges, forming n distinct bins
    edges = np.concatenate(([0.0], midpoints, [np.inf]))
    
    # 3. Digitize using these n+1 custom boundaries
    raw_indices = np.digitize(pressure_hpa, edges, right=True)
    
    # 4. Sort the dataset positions based on these bins
    sort_indices = np.argsort(raw_indices)
    binned_data = apply_filter(data, sort_indices)
    
    # 5. Sort mapped bin classifications
    sorted_raw_indices = raw_indices[sort_indices]
    
    # 6. Shift indices down by 1 so your n bins are 0-indexed (0 to n-1)
    bounded_indices = np.clip(sorted_raw_indices - 1, 0, len(bins) - 1)
    indices = bounded_indices

    # 7. Match dimensions: centers, labels, and heights will all be exactly size 18!
    centers = bins
    labels = bins  
    heights = np.diff(edges) 

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