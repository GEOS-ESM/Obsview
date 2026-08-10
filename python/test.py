import os
import re
import numpy as np
from netCDF4 import Dataset
from dataclasses import replace, dataclass
from datetime import datetime, timezone
from typing import Optional


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
    # time: Optional[np.ndarray]
    
    #error, also add other relevant variables used in calculations
    #scale
    kt: np.ndarray
    sid: np.ndarray
    #calculated variables
    amb: np.ndarray
    job: Optional[np.ndarray] = None
    joa: Optional[np.ndarray] = None
    esigo: Optional[np.ndarray] = None
    esigb: Optional[np.ndarray] = None
    #metadata: Metadata
    lev_type: Optional[str] = None #'pressure' or 'channel'
    all_lev: Optional[np.ndarray] = None
    fill_values: Optional[dict] = None
    datetime: Optional[object] = None
    

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
            "lon": nc.variables['lon'][:], 
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

            all_lev = np.unique(lev[lev< 1.0e15]),
        )
        return obj
    
    #Subfunction for turning array of seconds after epoch into single datetime object
    def _parse_datetime_from_filename(self, filename: str) -> datetime:
        base = os.path.basename(filename)

        # Match 'YYYYMMDD_HHz' (case-insensitive 'z').
        m = re.search(r"(\d{8})_(\d{2})z", base, flags=re.IGNORECASE)
        date_str, hour_str = m.group(1), m.group(2)

        # Build a UTC-aware datetime; strptime validates the calendar date.
        dt = datetime.strptime(date_str + hour_str, "%Y%m%d%H")
        return dt.replace(tzinfo=timezone.utc)    

    #Main reading method to be used to load and process ODS files
    def read(self, filename: str) -> ObservationData:
        nc = self._open_file(filename)
        raw = self._load_variables(nc)
        raw = self._calc_variables(raw)
        obj = self._flatten_data(raw)
        dt = self._parse_datetime_from_filename(filename)
        obj = replace(obj, datetime=dt)

        return obj



def main():
    reader = ODSReader()
    data = reader.read("python/x0050.diag_amsua_n19.20230801_00z.ods")
    print(f"Date time: {data.datetime}")


main()