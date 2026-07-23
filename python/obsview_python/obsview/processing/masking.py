#Module for creating masks of different data, each function returns a 1-D array of boolean values
import numpy as np
from ..loading.observationdata import ObservationData

#This should be for data that is valid, not necessarily passes qc
def valid_mask(data: ObservationData) -> np.ndarray: 
    missing_val = 1.0e15
    valid_mask = (
        (data.qc == 0)              #Should be changed later for missing value
        #& (data.sid == -999)        #Ditto
        & (data.kt == 40)           #This should be omitted later
        & (data.lev < missing_val)
        & (data.omb < missing_val)
        & (data.oma < missing_val)
        & (data.amb < missing_val)
    )
    return valid_mask

#This mask keeps data that isn't a missing value(same as valid_mask() but for data that contains specific fill values)
def fill_val_mask(data:ObservationData) -> np.ndarray:
    valid_mask = (
        (data.qc < np.abs(data.fill_values['qc']))                             
        & (data.omb < np.abs(data.fill_values['omb']))
        & (data.oma < np.abs(data.fill_values['oma']))
        & (data.sigo < np.abs(data.fill_values['sigo']))
        & (data.lev < np.abs(data.fill_values['lev']))

    )    
    return valid_mask

def qc_pass_mask(data: ObservationData) -> np.ndarray:
    qc_pass = (data.qc == 0)
    return qc_pass
   

def qc_fail_mask(data: ObservationData) -> np.ndarray:
    qc_fail = (data.qc > 0)
    return qc_fail