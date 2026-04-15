"""Configuration module for ODS handling."""
from typing import Any, Dict, List
from dataclasses import dataclass

@dataclass
class DataTypeInfo:
    value: int
    id: str
    units: str
    msigo: float

@dataclass
class DataSourceInfo:
    value: int
    id: str

@dataclass
class SensorInfo:
    id: str
    kxs: List[int]
    chns: List[int]

@dataclass
class QCFlagInfo:
    value: int
    id: str

class ODSConfig:
    def __init__(self):
        self._init_obs_attributes()
        self._init_data_types()
        self._init_data_sources()
        self._init_sensors()
        self._init_qc_flags()
    
    def _init_obs_attributes(self):
        self.obs_attributes = {'names': ['kx', 'ks', 'kt', 'time', 'lat', 'lon', 'lev', 'obs', 'sigo', 'omf', 'oma', 'qch', 'qcx', 'xm']}
    
    def _init_data_types(self):
        self.data_types = {44: DataTypeInfo(44, 'Upper-air virtual temperature', 'Kelvin', 10), 33: DataTypeInfo(33, 'Surface pressure', 'hPa', 12), 4: DataTypeInfo(4, 'Upper-air zonal wind', 'm/sec', 20), 5: DataTypeInfo(5, 'Upper-air meridional wind', 'm/sec', 20)}
        self.pressure_level_types = [4, 5, 11, 44]
        self.surface_types = [12, 33, 39]
        self.radiance_types = [40]
    
    def _init_data_sources(self):
        self.data_sources = {102: DataSourceInfo(102, 'SSM/I'), 120: DataSourceInfo(120, 'RAWINSONDE'), 220: DataSourceInfo(220, 'RAWINSONDE WINDS')}
    
    def _init_sensors(self):
        self.sensors = {'HIRS': SensorInfo('HIRS', [14, 16, 17], list(range(1, 20))), 'AMSUA': SensorInfo('AMSUA', [315, 316, 349], list(range(1, 16)))}
    
    def _init_qc_flags(self):
        self.qc_history_flags = {0: QCFlagInfo(0, 'none')}
        self.qc_exclusion_flags = {0: QCFlagInfo(0, 'none'), 1: QCFlagInfo(1, 'passive'), 2: QCFlagInfo(2, 'rejected by GSI')}
    
    def get_attribute(self, name):
        attrs = {'OBSATTRIBUTES': self.obs_attributes, 'KTS': self.data_types, 'KXS': self.data_sources, 'SENSORS': self.sensors, 'QCXS': self.qc_exclusion_flags}
        return attrs.get(name)

_config = ODSConfig()

def dconfig(*args):
    if not args:
        return {'OBSATTRIBUTES': _config.obs_attributes, 'KTS': _config.data_types, 'KXS': _config.data_sources, 'SENSORS': _config.sensors, 'QCXS': _config.qc_exclusion_flags}
    results = []
    for arg in args:
        value = _config.get_attribute(arg)
        if value is None:
            raise ValueError("Unknown configuration parameter: {}".format(arg))
        results.append(value)
    return results[0] if len(results) == 1 else tuple(results)
