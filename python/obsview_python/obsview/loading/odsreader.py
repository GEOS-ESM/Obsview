#Module containing ODSReader object 
import numpy as np
from netCDF4 import Dataset
from .observationdata import ObservationData


#TODO: add logic that populates lev_type with either "pressure" or "channel"


class ODSReader:

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
            "qc": nc.variables['qcexcl'][:],
            "lev": nc.variables['lev'][:],
            "kt": nc.variables['kt'][:],
            "sid": nc.variables['kx'][:],
            "lat": nc.variables['lat'][:],
            "lon": nc.variables['lon'][:]
        }
        return raw
    
     #Calculate new variables and append to raw dictionary
    def _calc_variables(self, raw: dict) -> dict:
        #Calculate
        amb = raw["omb"] - raw["oma"]
        #Append
        raw["amb"] = amb
        return raw

    #Flatten data, return ObservationData object
    def _flatten_data(self, raw: dict) -> ObservationData:
        lev = raw["lev"].flatten()

        obj = ObservationData(
            obs = raw["obs"].flatten(),
            omb = raw["omb"].flatten(),
            oma = raw["oma"].flatten(),
            sigo = raw["sigo"].flatten(),
            qc = raw["qc"].flatten(),
            lev = lev,
            lat = raw["lat"].flatten(),
            lon = raw['lon'].flatten(),

            kt = raw["kt"].flatten(),
            sid = raw["sid"].flatten(),

            amb = raw["amb"].flatten(),

            all_lev = np.unique(lev[lev< 1.0e15])
        )
        return obj

    #Main reading method to be used to load and process ODS files
    def read(self, filename: str) -> ObservationData:
        nc = self._open_file(filename)
        raw = self._load_variables(nc)
        raw = self._calc_variables(raw)
        obj = self._flatten_data(raw)

        return obj