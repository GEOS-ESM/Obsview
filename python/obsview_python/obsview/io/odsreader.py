#Module containing ODSReader object 
from netCDF4 import Dataset
from .observationdata import ObservationData

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
            "lev": nc.variables['lev'][:]
            #metadata to be added later
        }
        return raw
    
     #Calculate new variables, for now this is just for analysis minus background (amb)
     #and append to raw dictionary
    def _calc_variables(self, raw: dict) -> dict:
        amb = raw["omb"] - raw["oma"]
        raw["amb"] = amb
        return raw

    #Flatten data
    def _flatten_data(self, raw: dict) -> ObservationData:
        obj = ObservationData(
            obs = raw["obs"].flatten(),
            omb = raw["omb"].flatten(),
            oma = raw["oma"].flatten(),
            amb = raw["amb"].flatten(),
            sigo = raw["sigo"].flatten(),
            qc = raw["qc"].flatten(),
            lev = raw["lev"].flatten()
        )
        return obj

    #Main reading method to be used to load and process ODS files
    def read(self, filename: str) -> ObservationData:
        nc = self._open_file(filename)
        raw = self._load_variables(nc)
        raw = self._calc_variables(raw)
        obj = self._flatten_data(raw)

        return obj