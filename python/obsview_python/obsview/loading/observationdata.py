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
    obs: Optional[np.ndarray] = None
    omb: Optional[np.ndarray] = None
    omb_no_bias: Optional[np.ndarray] = None
    oma: Optional[np.ndarray] = None
    sigo: Optional[np.ndarray] = None
    qc: Optional[np.ndarray] = None
    bias: Optional[np.ndarray] = None
    lev: Optional[np.ndarray] = None

    lat: Optional[np.ndarray] = None
    lon: Optional[np.ndarray] = None
    #error
    #scale
    kt: Optional[np.ndarray] = None
    kx: Optional[np.ndarray]  = None     #sid for IODA files
    #calculated variables
    amb: Optional[np.ndarray] = None
    job: Optional[np.ndarray] = None
    joa: Optional[np.ndarray] = None
    esigo: Optional[np.ndarray] = None
    esigb: Optional[np.ndarray] = None
    #metadata: Metadata
    lev_type: Optional[str] = None          #'pressure' or 'channel'
    all_lev: Optional[np.ndarray] = None
    fill_values: Optional[dict] = None
    datetime: Optional[object] = None
    file_type: Optional[str] = None     #'ods' or 'ioda'
    


    #def concatenate()
    #def filter()
    #def group by time()
    #def average()
    #def subset()
    ...