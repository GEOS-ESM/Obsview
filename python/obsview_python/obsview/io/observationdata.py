#Module defining Metadata and a generic ObservationData classes

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
    amb: np.ndarray
    sigo: np.ndarray
    qc: np.ndarray
    lev: np.ndarray 
    #lev_type: pressure or channel
    #error, also add other relevant variables used in calculations

    #metadata: Metadata object

#For loading multiple files
@dataclass
class ObservationDataset:
    #def concatenate()
    #def filter()
    #def group by time()
    #def average()
    #def subset()
    ...