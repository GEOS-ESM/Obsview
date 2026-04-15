"""Data subsetting and filtering module."""
import numpy as np
from typing import Optional, Dict, Any, Union
from .models import ODSData

def odssubset(ods, criteria=None, **kwargs):
    if len(ods) == 0:
        return ods
    if isinstance(criteria, np.ndarray):
        return ods.subset(criteria)
    if isinstance(criteria, dict):
        mask = np.ones(len(ods), dtype=bool)
        for key, value in criteria.items():
            if hasattr(ods, key):
                attr = getattr(ods, key)
                if np.isscalar(value):
                    mask &= attr == value
                else:
                    mask &= np.isin(attr, value)
        return ods.subset(mask)
    if kwargs:
        mask = np.ones(len(ods), dtype=bool)
        if 'lat' in kwargs:
            lat_range = kwargs['lat']
            mask &= (ods.lat >= lat_range[0]) & (ods.lat <= lat_range[1])
        if 'lon' in kwargs:
            lon_range = kwargs['lon']
            mask &= (ods.lon >= lon_range[0]) & (ods.lon <= lon_range[1])
        if 'lev' in kwargs:
            lev_range = kwargs['lev']
            if np.isscalar(lev_range):
                mask &= ods.lev == lev_range
            else:
                mask &= (ods.lev >= lev_range[0]) & (ods.lev <= lev_range[1])
        return ods.subset(mask)
    return ods

def odsclean(ods):
    if len(ods) == 0:
        return ods
    mask = ods.qcx == 0
    return ods.subset(mask)
