#Module containing IODAReader class and relevant functions
import numpy as np
from netCDF4 import Dataset
from .observationdata import ObservationData

#TODO: add logic that populates lev_type with either "pressure" or "channel"

class IODAReader:
    
    #Open NetCDF file
    def _open_file(self,filename: str) -> Dataset:
        nc = Dataset(filename, "r")
        nc.set_auto_mask(False)
        return nc
    #Load data and flatten incoming arrays
    def _load_data(self, nc: Dataset) -> dict: 
        varname = "brightnessTemperature"     #Hardcoded for now, add function that takes user input to select variable name
        n_locations = np.size(nc.variables["Location"][:])
        raw = {
        "obs": nc.groups["ObsValue"].variables[varname][:].flatten(),
        "omb": nc.groups["ombg"].variables[varname][:].flatten(),
        "oma": nc.groups["oman"].variables[varname][:].flatten(),
        "sigo": nc.groups["EffectiveError0"].variables[varname][:].flatten(),
        "qc": nc.groups["EffectiveQC0"].variables[varname][:].flatten(),
        "all_lev": nc.variables["Channel"][:].flatten(),     #Hardcoded for now, change later to accept logic to determine what type of level variable(others include pressure and wavelength)
        "sid": 326,     #SID for Amsua Metop-B satellite, change later using config/rc file
        "kt": 40,       #Hardcoded for now, change later using config file
        "lev": np.tile(nc.variables["Channel"][:],n_locations)
        }
        return raw
        
    def _calc_variables(self, raw: dict) -> dict:
        #Calculate
        amb = raw["omb"] - raw["oma"]

        #Append
        raw["amb"] = amb
        return raw
    
    def _load_fill_values(self, nc: Dataset) -> dict:
        varname = "brightnessTemperature"  # keep consistent with _load_data

        # Map logical name -> the actual NetCDF variable object it was read from.
        # (Must mirror the sources used in _load_data.)
        var_sources = {
            "omb":  nc.groups["ombg"].variables[varname],
            "oma":  nc.groups["oman"].variables[varname],
            "sigo": nc.groups["EffectiveError0"].variables[varname],
            "qc":   nc.groups["EffectiveQC0"].variables[varname],
            "lev":  nc.variables["Channel"],
        }

        fill_values = {}
        for name, var in var_sources.items():
            if "_FillValue" in var.ncattrs():
                fill_values[name] = var.getncattr("_FillValue")
            else:
                fill_values[name] = None  # no declared fill value for this variable

        return fill_values   
    
    def _create_data_object(self, raw: dict, fill_values: dict) -> ObservationData:
        obj = ObservationData(
            obs = raw["obs"],
            omb = raw["omb"],
            oma = raw["oma"],
            sigo = raw["sigo"],
            qc = raw["qc"],
            lev = raw["lev"],
            kt = raw["kt"],
            sid = raw["sid"],
            amb = raw["amb"],
            all_lev= raw["all_lev"],
            fill_values = fill_values
        )
        return obj
        ...

    
    def read(self, filename: str) -> ObservationData:
        nc = self._open_file(filename)
        raw = self._load_data(nc)
        raw = self._calc_variables(raw)
        fill_values = self._load_fill_values(nc)
        obj = self._create_data_object(raw, fill_values)
        return obj