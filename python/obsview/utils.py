"""Utility functions for ODS handling."""
import numpy as np
from typing import Dict, Any
from .models import ODSData

def calculate_statistics(ods, attribute):
    data = getattr(ods, attribute)
    valid = np.isfinite(data)
    if not np.any(valid):
        return {'min': np.nan, 'max': np.nan, 'mean': np.nan, 'std': np.nan, 'median': np.nan, 'count': 0}
    return {'min': np.nanmin(data), 'max': np.nanmax(data), 'mean': np.nanmean(data), 'std': np.nanstd(data), 'median': np.nanmedian(data), 'count': np.sum(valid)}

def ods_to_dict(ods):
    return ods.to_dict()
