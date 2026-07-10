"""Small, generic helpers shared across the package."""
import numpy as np


def safe_filled(arr, fallback_value):
    """Fill masked values in an array.

    Floating-point arrays are filled with NaN (so downstream nan-aware
    statistics work correctly); everything else is filled with
    ``fallback_value``.
    """
    if np.issubdtype(arr.dtype, np.floating):
        return np.ma.filled(arr, fill_value=np.nan)
    else:
        return np.ma.filled(arr, fill_value=fallback_value)
