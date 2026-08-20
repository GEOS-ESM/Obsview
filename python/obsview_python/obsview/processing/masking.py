#Module for creating masks of different data, each function returns a 1-D array of boolean values
import numpy as np
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