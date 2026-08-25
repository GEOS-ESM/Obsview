#Module containing StatisticsData object
import numpy as np
from dataclasses import dataclass, field

@dataclass
class StatisticsData:
    nobs: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_omb: np.ndarray = field(default_factory = lambda: np.array([]))
    rms_omb: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_oma: np.ndarray = field(default_factory = lambda: np.array([]))
    rms_oma: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_job: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_joa: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_sigo: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_esigo: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_esigb: np.ndarray = field(default_factory = lambda: np.array([]))
    ...
