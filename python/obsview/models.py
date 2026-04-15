"""Data models for ODS observation data."""
import numpy as np
from typing import Optional, List, Dict, Any
from dataclasses import dataclass, field

@dataclass
class ODSData:
    kt: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int8))
    kx: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int16))
    ks: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int32))
    lon: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))
    lat: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))
    lev: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))
    time: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))
    obs: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))
    omf: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))
    oma: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))
    xm: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))
    sigo: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.float32))
    qch: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int16))
    qcx: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int8))
    filename: Optional[str] = None
    first_julian_day: int = 0
    latest_julian_day: int = 0
    latest_synoptic_hour: int = 0
    synoptic_hours_per_day: int = 0
    version: str = ""
    kt_names: Optional[np.ndarray] = None
    kt_units: Optional[np.ndarray] = None
    kx_names: Optional[np.ndarray] = None
    kx_meta: Optional[np.ndarray] = None
    qcx_names: Optional[np.ndarray] = None
    _scales: Dict[str, float] = field(default_factory=dict)
    _offsets: Dict[str, float] = field(default_factory=dict)
    _missing_values: Dict[str, float] = field(default_factory=dict)
    cidx: Optional[np.ndarray] = None
    cinfo: Optional[List] = None
    ssid: Optional[str] = None
    sslat: Optional[np.ndarray] = None
    sslon: Optional[np.ndarray] = None
    sslev: Optional[np.ndarray] = None
    sskt: Optional[np.ndarray] = None
    sskx: Optional[np.ndarray] = None
    ssqcx: Optional[np.ndarray] = None
    ssqch: Optional[np.ndarray] = None
    
    def __len__(self):
        return len(self.kt) if self.kt is not None else 0
    
    def __repr__(self):
        fname = Path(self.filename).name if self.filename else "None"
        return "ODSData(nobs={}, file={})".format(len(self), fname)
    
    def subset(self, mask):
        subset = ODSData()
        for attr in ['kt', 'kx', 'ks', 'lon', 'lat', 'lev', 'time', 'obs', 'omf', 'oma', 'xm', 'sigo', 'qch', 'qcx']:
            arr = getattr(self, attr)
            if arr is not None and len(arr) > 0:
                setattr(subset, attr, arr[mask])
        subset.filename = self.filename
        subset.first_julian_day = self.first_julian_day
        subset.latest_julian_day = self.latest_julian_day
        subset._scales = self._scales.copy()
        subset._offsets = self._offsets.copy()
        subset._missing_values = self._missing_values.copy()
        return subset
    
    def to_dict(self):
        return {'kt': self.kt, 'kx': self.kx, 'ks': self.ks, 'lon': self.lon, 'lat': self.lat, 'lev': self.lev, 'time': self.time, 'obs': self.obs, 'omf': self.omf, 'oma': self.oma, 'xm': self.xm, 'sigo': self.sigo, 'qch': self.qch, 'qcx': self.qcx}
