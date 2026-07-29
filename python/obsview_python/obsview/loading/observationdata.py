#Module defining Metadata and a generic ObservationData classes
from typing import Optional
from dataclasses import dataclass
import numpy as np
@dataclass
class Metadata:
    kt: np.ndarray
    kx: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    datetime: np.ndarray
    filename: np.ndarray
    version: np.ndarray

@dataclass
class ObservationData:
    obs: np.ndarray
    omb: np.ndarray
    oma: np.ndarray
    sigo: np.ndarray
    qc: np.ndarray
    lev: np.ndarray

    lat: Optional[np.ndarray] 
    lon: Optional[np.ndarray]
    # time: Optional[np.ndarray]
    
    #error, also add other relevant variables used in calculations
    #scale
    kt: np.ndarray
    sid: np.ndarray
    #calculated variables
    amb: np.ndarray
    job: Optional[np.ndarray] = None
    joa: Optional[np.ndarray] = None
    esigo: Optional[np.ndarray] = None
    esigb: Optional[np.ndarray] = None
    #metadata: Metadata
    lev_type: Optional[str] = None #'pressure' or 'channel'
    all_lev: Optional[np.ndarray] = None
    fill_values: Optional[dict] = None
    datetime: Optional[object] = None
    
#For loading multiple files
@dataclass
class ObservationDataset:
    #def concatenate()
    #def filter()
    #def group by time()
    #def average()
    #def subset()
    ...