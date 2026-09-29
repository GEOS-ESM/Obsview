#Miscellaneous functions
import time
import numpy as np
from datetime import datetime, timezone


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

#Timer function for performance benchmarking
def timer(base_fn):
    def enhanced_fn():
        start_time = time.perf_counter()
        base_fn()
        end_time = time.perf_counter()
        print(f"Task time: {end_time - start_time} seconds")
    return enhanced_fn  

def _parse_time(s: str) -> datetime:
    # Accept 'YYYY-MM-DD HH
    dt = datetime.strptime(s, "%Y-%m-%d %H")
    return dt.replace(tzinfo=timezone.utc)
