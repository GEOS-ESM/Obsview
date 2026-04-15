"""ODS data loading module."""
import numpy as np
import xarray as xr
from typing import Optional, Union, List
from pathlib import Path
import logging
from .models import ODSData

logger = logging.getLogger(__name__)

def odsload(filename, jdays=None, hours=None, attrs=None):
    file_path = Path(filename)
    if not file_path.exists():
        raise FileNotFoundError("ODS file not found: {}".format(filename))
    logger.info("Loading ODS file: {}".format(filename))
    ods = _get_ods_info(filename)
    if attrs is None:
        attrs = ['kt', 'kx', 'ks', 'lon', 'lat', 'lev', 'time', 'obs', 'omf', 'oma', 'xm', 'qcx', 'qch', 'sigo']
    logger.info("Successfully loaded {} observations".format(len(ods)))
    return ods

def _get_ods_info(odsfile):
    file_path = Path(odsfile)
    if not file_path.exists():
        raise FileNotFoundError("ODS file not found: {}".format(odsfile))
    try:
        ds = xr.open_dataset(odsfile)
    except Exception as e:
        raise ValueError("Cannot open file as NetCDF: {}".format(odsfile))
    try:
        ods = ODSData()
        ods.filename = odsfile
        for attr in ['first_julian_day', 'latest_julian_day', 'latest_synoptic_hour', 'version']:
            if attr in ds.attrs:
                setattr(ods, attr, ds.attrs[attr])
        if not ods.version or float(str(ods.version)[0]) < 2:
            raise ValueError("{}: Not an ODS Version 2 file.".format(odsfile))
        if 'syn_beg' in ds:
            ods.synoptic_hours_per_day = ds['syn_beg'].dims[0] if ds['syn_beg'].dims else 4
        logger.info("Read ODS header: JD {}-{}".format(ods.first_julian_day, ods.latest_julian_day))
        return ods
    finally:
        ds.close()
