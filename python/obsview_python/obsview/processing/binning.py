#Module containing BinnedData object as well as functions to create bins and sort data by level
import numpy as np
from typing import Optional, List
from dataclasses import dataclass
from ..loading.observationdata import ObservationData
from ..reference import Reference, Observation
from .filtering import apply_filter

@dataclass
class BinnedData:
    data: ObservationData           #Rearranged data by bin
    bin_indices: Optional[np.ndarray] = None   
    bin_labels: Optional[np.ndarray] = None   
    #This should eventually replace the large ObservationData attribute
    nobs: Optional[np.ndarray] = None       #array with number of observations (int) for each vertical level 
    is_ts: Optional[bool] = None         #For plotting, tells us if the binned data came from a time series object
    ts_range: Optional[List[object]] = None   #For plotting, gives start and end datetimes










def create_bins(data: ObservationData, obs: Observation, ref: Reference) -> BinnedData:
    if data.lev_type == "channel":
        return create_channel_bins(data)
    elif data.lev_type == "pressure":
        #Get pressure level set    
        bins = np.array(ref.get_pressure_levels(obs))
        return create_pressure_bins(data, bins)
    raise ValueError(f"Unknown lev_type: {data.lev_type!r}")





def create_pressure_bins(data: ObservationData, bins = np.ndarray) -> BinnedData:
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
    
    # 6. Shift indices down by 1 so n bins are 0-indexed (0 to n-1)
    bounded_indices = np.clip(sorted_raw_indices - 1, 0, len(bins) - 1)
    indices = bounded_indices

    labels = bins  

    obj = BinnedData(
        data=binned_data,
        bin_indices=indices,
        bin_labels=labels
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
         

    obj = BinnedData(
        data=binned_data,
        bin_indices=indices,
        bin_labels=channels
    )
    return obj