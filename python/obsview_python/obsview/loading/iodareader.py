#Module containing IODAReader class and relevant functions
import numpy as np
from netCDF4 import Dataset
from dataclasses import replace
from datetime import datetime, timezone
from .observationdata import ObservationData
from ..config import VARNAME_TO_KT


class IODAReader:
    
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
    #Load data and flatten incoming arrays
    def _load_data(self, nc: Dataset, varname: str, kx: int) -> dict: 
        lev_type = self._get_lev_type(varname)
        n_locations = np.size(nc.variables["Location"][:])

        raw = {
        "obs": nc.groups["ObsValue"].variables[varname][:].flatten(),
        "omb": nc.groups["ombg"].variables[varname][:].flatten(),
        "oma": nc.groups["oman"].variables[varname][:].flatten(),
        "sigo": nc.groups["EffectiveError0"].variables[varname][:].flatten(),
        "qc": nc.groups["EffectiveQC0"].variables[varname][:].flatten(),
        "datetime": nc.groups["MetaData"].variables["dateTime"][:].flatten(),
        "bias": nc.groups["ObsBias1"].variables[varname][:].flatten(),
        
        "kx": kx,     #No better way to retrieve kx/sid for now
        "kt": VARNAME_TO_KT.get(varname)       
        }

    #Level-type specific variables 
        if lev_type == 'pressure':
            raw["lev"] = nc.groups["MetaData"].variables["pressure"][:]
            raw["all_lev"] = np.unique(raw["lev"])
            raw["lat"] = nc.groups["MetaData"].variables["latitude"][:]
            raw["lon"] = nc.groups["MetaData"].variables["longitude"][:]
            ...
        elif lev_type == 'channel':
            n_channels = np.size(nc.variables["Channel"][:])
            raw["all_lev"] = nc.variables["Channel"][:].flatten()
            raw["lev"] = np.tile(nc.variables["Channel"][:],n_locations)
            raw["lat"] = np.repeat(nc.groups["MetaData"].variables["latitude"][:], n_channels)
            raw["lon"] = np.repeat(nc.groups["MetaData"].variables["longitude"][:], n_channels)
        else:
            raise ValueError(f"Unknown lev_type: {lev_type!r}")


        return raw
        

        
    def _calc_variables(self, raw: dict) -> dict:
        #Calculate
        amb = raw["omb"] - raw["oma"]
        omb_no_bias = raw["omb"]+raw["bias"]
        #Append
        raw["amb"] = amb
        raw["omb_no_bias"] = omb_no_bias
        return raw
    
    def _load_fill_values(self, nc: Dataset, varname: str) -> dict:
    
        var_sources = {
            "omb":  nc.groups["ombg"].variables[varname],
            "oma":  nc.groups["oman"].variables[varname],
            "sigo": nc.groups["EffectiveError0"].variables[varname],
            "qc":   nc.groups["EffectiveQC0"].variables[varname],
            #"lev":  nc.variables["Channel"],
            "lat": nc.groups['MetaData'].variables['latitude'],
            "lon": nc.groups['MetaData'].variables['longitude']
        }

        fill_values = {}
        for name, var in var_sources.items():
            if "_FillValue" in var.ncattrs():
                fill_values[name] = var.getncattr("_FillValue")
            else:
                fill_values[name] = None  # no declared fill value for this variable

        return fill_values 

    #Subfunction for turning array of seconds after epoch into single datetime object
    def _synoptic_time_from_datetimes(self, epoch_seconds: np.ndarray) -> datetime:
   
    # Guard against fill values / non-finite entries before taking the median.
        fill = -9223372036854775801
        valid = epoch_seconds[np.isfinite(epoch_seconds)]
        if fill is not None:
            valid = valid[valid != fill]
        if valid.size == 0:
            raise ValueError("No valid dateTime values to determine synoptic time.")

        # Median epoch -> center of the observation window.
        median_epoch = float(np.median(valid))

        # Round to the nearest 6-hour boundary (6h = 21600 s).
        six_hours = 6 * 3600
        rounded_epoch = round(median_epoch / six_hours) * six_hours

        # Build a timezone-aware UTC datetime.
        return datetime.fromtimestamp(rounded_epoch, tz=timezone.utc)  
    
    def _create_data_object(self, raw: dict, fill_values: dict, varname: str) -> ObservationData:
        level_type = self._get_lev_type(varname)
        obj = ObservationData(
            obs = raw["obs"],
            omb = raw["omb"],
            omb_no_bias = raw["omb_no_bias"],
            oma = raw["oma"],
            sigo = raw["sigo"],
            qc = raw["qc"],
            lev = raw["lev"],
            lat = raw["lat"],
            lon = raw["lon"],
            kt = raw["kt"],
            kx = raw["kx"],
            amb = raw["amb"],
            all_lev= raw["all_lev"],
            fill_values = fill_values,
            lev_type = level_type,
            file_type = 'ioda' 
        )
        return obj
        ...

    
    def read(self, filename: str, varname: str, kx: int) -> ObservationData:

        nc = self._open_file(filename)
        raw = self._load_data(nc, varname, kx)
        raw = self._calc_variables(raw)
        fill_values = self._load_fill_values(nc, varname)
        obj = self._create_data_object(raw, fill_values, varname)
        #Make datetime object
        synoptic = self._synoptic_time_from_datetimes(raw["datetime"])
        obj = replace(obj, datetime=synoptic)
        return obj

