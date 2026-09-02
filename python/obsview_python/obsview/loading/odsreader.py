#Module containing ODSReader object 
import os
import re
import numpy as np
from netCDF4 import Dataset
from dataclasses import replace
from datetime import datetime, timezone
from .observationdata import ObservationData
from ..config import VARNAME_TO_KT
from ..processing.filtering import apply_filter


class ODSReader:

    def _get_lev_type(self, varname: str) -> str:
        if varname == 'brightnessTemperature':
            lev_type = 'channel'
        else:
            lev_type = 'pressure'
        return lev_type
    
    #Open NetCDF file
    def _open_file(self,filename: str) -> Dataset:
        nc = Dataset(filename, "r")
        nc.set_auto_mask(False)
        return nc

    #Load variable data into Observation data class
    def _load_variables(self, nc: Dataset) -> dict:
        raw = {
            "obs": nc.variables['obs'][:],
            "omb": nc.variables['omf'][:],
            "oma": nc.variables['oma'][:],
            "sigo": nc.variables['xvec'][:],
            "bias": nc.variables['xm'][:],
            "qc": nc.variables['qcexcl'][:],
            "lev": nc.variables['lev'][:],
            "kt": nc.variables['kt'][:],
            "kx": nc.variables['kx'][:],
            "lat": nc.variables['lat'][:],
            "lon": nc.variables['lon'][:], 
        }
        return raw
    
     #Calculate new variables and append to raw dictionary
    def _calc_variables(self, raw: dict) -> dict:
        #Calculate
        amb = raw["omb"] - raw["oma"]
        omb_no_bias = raw["omb"] + raw['bias']
        #Append
        raw["amb"] = amb
        raw["omb_no_bias"] = omb_no_bias
        return raw

    #Flatten data, return ObservationData object
    def _flatten_data(self, raw: dict, varname: str) -> ObservationData:
        lev = raw["lev"].flatten()
        level_type = self._get_lev_type(varname)
        obj = ObservationData(
            obs = raw["obs"].flatten(),
            omb = raw["omb"].flatten(),
            oma = raw["oma"].flatten(),
            sigo = raw["sigo"].flatten(),
            qc = raw["qc"].flatten(),
            lev = lev,
            lat = raw["lat"].flatten(),
            lon = raw['lon'].flatten(),
            bias = raw["bias"].flatten(),

            kt = raw["kt"].flatten(),
            kx = raw["kx"].flatten(),

            amb = raw["amb"].flatten(),
            omb_no_bias = raw["omb_no_bias"].flatten(),

            all_lev = np.unique(lev[lev< 1.0e15]),
            lev_type = level_type,
            file_type = 'ods'
        )
        return obj
    
    #Create a mask to keep data with a unique kt and kx
    def _kt_kx_mask(self, obj: ObservationData, varname: str, kx: int) -> np.ndarray:
        kt = VARNAME_TO_KT.get(varname)
        valid_mask = ((obj.kt == kt)
                      & (obj.kx == kx))
        return valid_mask
    
    def _filter_kt_and_kx(self, obj: ObservationData, mask: np.ndarray) -> ObservationData:
        obj = apply_filter(obj, mask)
        return obj

    #Subfunction for parsing filename to provide a datetime object
    def _parse_datetime_from_filename(self, filename: str) -> datetime:
        base = os.path.basename(filename)

        # Match 'YYYYMMDD_HHz' (case-insensitive 'z').
        m = re.search(r"(\d{8})_(\d{2})z", base, flags=re.IGNORECASE)
        date_str, hour_str = m.group(1), m.group(2)

        # Build a UTC-aware datetime; strptime validates the calendar date.
        dt = datetime.strptime(date_str + hour_str, "%Y%m%d%H")
        return dt.replace(tzinfo=timezone.utc)    

    #Main reading method to be used to load and process ODS files
    def read(self, filename: str, varname: str, kx: int) -> ObservationData:
        nc = self._open_file(filename)
        raw = self._load_variables(nc)
        raw = self._calc_variables(raw)
        obj = self._flatten_data(raw, varname)
        kt_mask = self._kt_kx_mask(obj, varname, kx)
        obj = self._filter_kt_and_kx(obj, kt_mask)
        dt = self._parse_datetime_from_filename(filename)
        single_kt = obj.kt[0]
        single_kx = obj.kx[0]
        obj = replace(obj, datetime=dt, kx = single_kx, kt = single_kt)
        
        

        return obj
