#Module for creating masks of different data, each function returns a 1-D array of boolean values
import numpy as np
from typing import Tuple
from ..loading.observationdata import ObservationData


#This mask keeps data that isn't a missing value(same as valid_mask() but for data that contains specific fill values)
def fill_val_mask(data:ObservationData) -> np.ndarray:
    if data.fill_values == None:
        qc_fill_val = -127
        fill_val = 1.0e15
        valid_mask = (
            (data.qc != qc_fill_val)
            &(data.lev < fill_val)
        )
    else:
        valid_mask = (
            (data.qc < np.abs(data.fill_values['qc']))                             
            & (data.omb < np.abs(data.fill_values['omb']))
            & (data.oma < np.abs(data.fill_values['oma']))
            & (data.sigo < np.abs(data.fill_values['sigo']))
        )    
    return valid_mask


#Latitude/longitude masks

def valid_latlon_mask(data: ObservationData) -> np.ndarray:
    valid_mask = (
        (data.lat < np.abs(data.fill_values["lat"]))
        &(data.lon < np.abs(data.fill_values["lat"]))
    )
    return valid_mask


# def nh_mask
# def sh_mask
# def tr_mask





#Quality control masks
def qc_pass_mask(data: ObservationData) -> np.ndarray:
    qc_pass = (data.qc == 0)
    return qc_pass
   

def qc_fail_mask(data: ObservationData) -> np.ndarray:
    qc_fail = (data.qc > 0)
    return qc_fail

#Masks of u and v wind components to match observations
def uv_match_masks(u: ObservationData, v: ObservationData) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return index arrays (iu, iv) such that u[iu] and v[iv] are the same
    observation.
    """
    #Fast way (two sets are already aligned)
    if (u.lat.shape == v.lat.shape
            and np.array_equal(u.lat, v.lat)
            and np.array_equal(u.lon, v.lon)
            and np.array_equal(u.lev, v.lev)):
        idx = np.arange(u.lat.size)
        return idx, idx
    
    ku = _keys(u)
    kv = _keys(v)
    _, iu, iv = np.intersect1d(ku, kv, return_indices=True)

    # dropped = u.lat.size - iu.size
    # if dropped:
    #     log.warning("windVector: %d of %d U records had no matching V record "
    #                 "(or were duplicate keys)", dropped, u.lat.size)
    return iu, iv


def _keys(d: ObservationData) -> np.ndarray:
    decimals = 3    #Number of decimals to round observations


    rows = np.column_stack([np.round(d.lat, decimals),
                            np.round(d.lon, decimals),
                            np.round(d.lev, decimals)]).astype(np.float64)
    rows = np.ascontiguousarray(rows)
    return rows.view(np.dtype((np.void, rows.dtype.itemsize * 3))).ravel()

    ...