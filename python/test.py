import os
import re
import numpy as np
from netCDF4 import Dataset
from dataclasses import replace, dataclass
from datetime import datetime, timezone
from typing import Optional


#config.py
VARNAME_TO_KT = {
    "bendingAngle": 89,
    "windEastward": 4,
    "windNorthward": 5,
    "specificHumidity": 11,
    "virtualTemperature": 44,
    "airTemperature": 44,
    "ozoneProfile": 87,
    "brightnessTemperature": 40,
}







@dataclass
class ObservationData:
    obs: np.ndarray
    omb: np.ndarray
    oma: np.ndarray
    sigo: np.ndarray
    qc: np.ndarray
    lev: np.ndarray

    lat: Optional[np.ndarray] 
    lon: Optional[np.ndarray]
    #error
    #scale
    kt: Optional[np.ndarray]
    sid: Optional[np.ndarray]
    #calculated variables
    amb: Optional[np.ndarray]
    job: Optional[np.ndarray] = None
    joa: Optional[np.ndarray] = None
    esigo: Optional[np.ndarray] = None
    esigb: Optional[np.ndarray] = None
    #metadata: Metadata
    lev_type: Optional[str] = None          #'pressure' or 'channel'
    all_lev: Optional[np.ndarray] = None
    fill_values: Optional[dict] = None
    datetime: Optional[object] = None
    



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
    def _load_data(self, nc: Dataset, varname: str) -> dict: 
        lev_type = self._get_lev_type(varname)
        n_locations = np.size(nc.variables["Location"][:])

        raw = {
        "obs": nc.groups["ObsValue"].variables[varname][:].flatten(),
        "omb": nc.groups["ombg"].variables[varname][:].flatten(),
        "oma": nc.groups["oman"].variables[varname][:].flatten(),
        "sigo": nc.groups["EffectiveError0"].variables[varname][:].flatten(),
        "qc": nc.groups["EffectiveQC0"].variables[varname][:].flatten(),
        "datetime": nc.groups["MetaData"].variables["dateTime"][:].flatten(),

        #TODO: Change these to not be hardcoded
        "sid": nc.groups["MetaData"].variables["satelliteIdentifier"][:].flatten(),
        "kt": 40       
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
        #Append
        raw["amb"] = amb
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
            oma = raw["oma"],
            sigo = raw["sigo"],
            qc = raw["qc"],
            lev = raw["lev"],
            lat = raw["lat"],
            lon = raw["lon"],
            kt = raw["kt"],
            sid = raw["sid"],
            amb = raw["amb"],
            all_lev= raw["all_lev"],
            fill_values = fill_values,
            lev_type = level_type 
        )
        return obj
        ...

    
    def read(self, filename: str, varname: str) -> ObservationData:

        nc = self._open_file(filename)
        raw = self._load_data(nc, varname)
        raw = self._calc_variables(raw)
        fill_values = self._load_fill_values(nc, varname)
        obj = self._create_data_object(raw, fill_values, varname)
        #Make datetime object
        synoptic = self._synoptic_time_from_datetimes(raw["datetime"])
        obj = replace(obj, datetime=synoptic)
        return obj


#binning.py



filename = "python/data/IODA files/j54rp1.jedi_hofx.20260101_03z/sondes.20260101T030000Z.nc4"
#filename = "python/amsua_metop-b.20260125T150000Z.nc4"
varname = "windEastward"        #(U) zonal wind
#varname = "brightnessTemperature"


def main():
    reader = IODAReader()
    data = reader.read(filename, varname)
    print(f"Unique sid: {np.unique(data.sid)}")

    


main()