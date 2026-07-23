import time
from typing import Dict, Optional
from dataclasses import dataclass, fields, field, replace

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt


###################################
#         Loading Package         #
###################################

#observationdata.py
@dataclass
class Metadata:
    #kt: np.ndarray
    # kx: np.ndarray
    # lat: np.ndarray
    # lon: np.ndarray
    datetime: Optional[np.ndarray]
    filename: Optional[np.ndarray]
    version: Optional[np.ndarray]
    kt_names: Optional[np.ndarray]
    kx_names: Optional[np.ndarray]

@dataclass
class ObservationData:
    obs: np.ndarray
    omb: np.ndarray
    oma: np.ndarray
    sigo: np.ndarray
    qc: np.ndarray
    lev: np.ndarray

    # lat: Optional[np.ndarray] #optional for now
    # lon: Optional[np.ndarray]
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

#odsreader
@dataclass
class ObservationDataset:
    #def concatenate()
    #def filter()
    #def group by time()
    #def average()
    #def subset()
    ...

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
            "sid": nc.variables['kx'][:]
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

            kt = raw["kt"].flatten(),
            sid = raw["sid"].flatten(),

            amb = raw["amb"].flatten(),
            job = raw["job"].flatten(),
            joa = raw["joa"].flatten(),
            esigo = raw["esigo"].flatten(),
            esigb = raw["esigb"].flatten(),

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
#TODO: add logic that populates lev_type with either "pressure" or "channel"

#iodareader.py
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
#TODO: add logic that populates lev_type with either "pressure" or "channel"








###############################
#       Processing Package    #
###############################

#masking.py

#This should be for data that is valid, not necessarily passes qc
def valid_mask(data: ObservationData) -> np.ndarray: 
    missing_val = 1.0e15
    valid_mask = (
        (data.qc == 0)              #Should be changed later for missing value
        #& (data.sid == -999)        #Ditto
        & (data.kt == 40)           #This should be omitted later
        & (data.lev < missing_val)
        & (data.omb < missing_val)
        & (data.oma < missing_val)
        & (data.amb < missing_val)
    )
    return valid_mask

#This mask keeps data that isn't a missing value(same as valid_mask() but for data that contains specific fill values)
def fill_val_mask(data:ObservationData) -> np.ndarray:
    valid_mask = (
        (data.qc < np.abs(data.fill_values['qc']))                             
        & (data.omb < np.abs(data.fill_values['omb']))
        & (data.oma < np.abs(data.fill_values['oma']))
        & (data.sigo < np.abs(data.fill_values['sigo']))
        & (data.lev < np.abs(data.fill_values['lev']))

    )    
    return valid_mask

def qc_pass_mask(data: ObservationData) -> np.ndarray:
    qc_pass = (data.qc == 0)
    return qc_pass
   

def qc_fail_mask(data: ObservationData) -> np.ndarray:
    qc_fail = (data.qc > 0)
    return qc_fail
    



#filtering.py
def apply_filter(data: ObservationData, mask: np.ndarray) -> ObservationData:
    updated_fields = {}
    if np.shape(data.obs) != np.shape(data.lev):
        ignore_fields = {"lev","all_lev"}
    else:
        ignore_fields = {"all_lev"}

    for field in fields(data):
        value = getattr(data, field.name)

        if isinstance(value, np.ndarray) and field.name not in ignore_fields: #If the field is type ndarray and not in the ignore fields dict
            updated_fields[field.name] = value[mask]
        else:
            updated_fields[field.name] = value

    return replace(data, **updated_fields)


#derived.py  Module for calculating derived values like job, joa, esigo, and esigb
#This prevents overflow errors from trying to calculate these values before masking
#since fill values are large numbers(~1e38)
def calc_derived(data: ObservationData) -> ObservationData:
    """
    Compute job, joa, esigo, esigb (and amb if missing) from omb/oma/sigo.
    Must be called AFTER filtering so no fill values remain (prevents overflow).
    """
    amb = data.amb if data.amb is not None else (data.omb - data.oma)

    job = data.omb**2 / data.sigo**2
    joa = data.oma**2 / data.sigo**2
    esigo = data.omb * data.oma
    esigb = data.omb * amb

    return replace(
        data,
        amb=amb,
        job=job,
        joa=joa,
        esigo=esigo,
        esigb=esigb,
    )


#binning.py
@dataclass
class BinnedData:
    data: ObservationData           #Rearranged data by bin
    bin_centers: np.ndarray         #Averaged value between bin levels (used for plotting)
    bin_indices: np.ndarray         
    bin_labels: np.ndarray          #Array of each bin level (unique)
    bin_heights: np.ndarray
    #level_type: str (pressure or channel)

#TODO: define this function    
def create_pressure_bins():
    ...

def create_channel_bins(data: ObservationData) -> BinnedData:
    channels = np.unique(data.all_lev)   # Same as bin_labels
    
    # Create a sorting index that arranges data by channel
    sort_indices = np.argsort(data.lev)
    
    # Apply the sorting to get data organized by bin
    binned_data = apply_filter(data, sort_indices)
    
    # Now calculate bin indices based on the sorted data
    indices = np.searchsorted(channels, binned_data.lev)
    
    n = len(channels)
    centers = np.arange(1, n + 1)            
    heights = 0.95 * np.ones(n)      

    obj = BinnedData(
        data=binned_data,
        bin_centers=centers,
        bin_indices=indices,
        bin_labels=channels,
        bin_heights=heights
    )
    return obj




#####################################################
#                   Stats                           #
#####################################################

#statisticsdata.py
@dataclass
class StatisticsData:
    nobs: np.ndarray = field(default_factory = lambda: np.array([]))
    #nonobs, to be used later for plotting unused observations (red bars)
    mean_omb: np.ndarray = field(default_factory = lambda: np.array([]))
    rms_omb: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_oma: np.ndarray = field(default_factory = lambda: np.array([]))
    rms_oma: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_job: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_joa: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_sigo: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_esigo: np.ndarray = field(default_factory = lambda: np.array([]))
    mean_esigb: np.ndarray = field(default_factory = lambda: np.array([]))
    ...

#calculate_stats.py
def count_obs_per_bin(binned_data: BinnedData) -> np.ndarray:
    """Return per-bin observation counts only (no derived-variable access)."""
    n_bins = len(binned_data.bin_labels)
    return np.bincount(binned_data.bin_indices, minlength=n_bins)



def calculate_stats(binned_data: BinnedData) -> StatisticsData:
    n_bins = len(binned_data.bin_labels)
    
    # Count observations per bin
    nobs = np.bincount(binned_data.bin_indices, minlength=n_bins)
    
    # Find the start index of each bin in the sorted data
    bin_starts = np.searchsorted(binned_data.bin_indices, np.arange(n_bins))
    bin_ends = np.append(bin_starts[1:], len(binned_data.bin_indices))
    
    # Pre-allocate arrays
    mean_omb = np.full(n_bins, np.nan)
    mean_oma = np.full(n_bins, np.nan)
    rms_omb = np.full(n_bins, np.nan)
    rms_oma = np.full(n_bins, np.nan)
    mean_job = np.full(n_bins, np.nan)
    mean_joa = np.full(n_bins, np.nan)
    mean_sigo = np.full(n_bins, np.nan)
    mean_esigo = np.full(n_bins, np.nan)
    mean_esigb = np.full(n_bins, np.nan)
    
    # Calculate statistics for each bin
    for i in range(n_bins):
        if nobs[i] > 0:
            start, end = bin_starts[i], bin_ends[i]
            
            mean_omb[i] = np.mean(binned_data.data.omb[start:end])
            mean_oma[i] = np.mean(binned_data.data.oma[start:end])
            rms_omb[i] = np.sqrt(np.mean(binned_data.data.omb[start:end]**2))
            rms_oma[i] = np.sqrt(np.mean(binned_data.data.oma[start:end]**2))
            mean_job[i] = np.mean(binned_data.data.job[start:end])
            mean_joa[i] = np.mean(binned_data.data.joa[start:end])
            mean_sigo[i] = np.mean(binned_data.data.sigo[start:end])
            mean_esigo[i] = np.sqrt(np.abs(np.mean(binned_data.data.esigo[start:end])))
            mean_esigb[i] = np.sqrt(np.abs(np.mean(binned_data.data.esigb[start:end])))
    
    return StatisticsData(
        nobs=nobs,
        mean_omb=mean_omb,
        mean_oma=mean_oma,
        rms_omb=rms_omb,
        rms_oma=rms_oma,
        mean_job=mean_job,
        mean_joa=mean_joa,
        mean_sigo=mean_sigo,
        mean_esigo=mean_esigo,
        mean_esigb=mean_esigb
    )
    ...


######################################
#           Plotting Package         #
######################################
#panels.py
#TODO: Change y labels to reflect lev_type
def plot_stats(pass_data: BinnedData, fail_data:BinnedData, stats: StatisticsData):
    
    fig = plt.figure(figsize=(10, 7))  # hardcoded size for now

    plt.subplot(2, 2, 1)
    _panel_nobs(pass_data, fail_data)

    plt.subplot(2, 2, 2)
    _panel_resstats(pass_data, stats)

    plt.subplot(2, 2, 3)
    _panel_jo(pass_data, stats)

    plt.subplot(2, 2, 4)
    _panel_sigo(pass_data, stats)

    return fig


def _panel_nobs(pass_data: BinnedData, fail_data: BinnedData):
    
    pass_nobs = count_obs_per_bin(pass_data)
    fail_nobs = count_obs_per_bin(fail_data)

    bin_centers = pass_data.bin_centers
    bin_heights = pass_data.bin_heights
    bar_width = bin_heights * 0.8   # wider single bar since we overlap now

    for i in range(len(bin_centers)):
        y = bin_centers[i]

        # Draw the "not used" (fail) bar first, fully opaque.
        plt.barh(
            y, fail_nobs[i],
            height=bar_width[i],
            color="red",
            label="not used" if i == 0 else "",
            zorder=1,
        )

        # Draw the "used" (pass) bar on top, at the SAME y, with opacity
        # so the red underneath is still visible.
        plt.barh(
            y, pass_nobs[i],
            height=bar_width[i],
            color="green",
            alpha=0.6,
            label="used" if i == 0 else "",
            zorder=2,
        )

    plt.ylabel("Channel")
    plt.title("Observation Count vs Channel")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()


def _panel_resstats(binned_data: BinnedData, stats: StatisticsData):
    """Panel 2: Mean & RMS of Obs Residuals vs Channel."""
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels = ["Mean o-b", "Mean o-a", "RMS o-b", "RMS o-a"]
    labeled = {"omb": False, "oma": False, "rms_omb": False, "rms_oma": False}

    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], stats.mean_omb[i], height=bar_width[i],
                 color="cyan", label=labels[0] if not labeled["omb"] else "")
        labeled["omb"] = True

        plt.barh(y + offsets[2] * bar_width[i], stats.mean_oma[i], height=bar_width[i],
                 color="orange", label=labels[1] if not labeled["oma"] else "")
        labeled["oma"] = True

        plt.barh(y + offsets[1] * bar_width[i], stats.rms_omb[i], height=bar_width[i],
                 color="blue", label=labels[2] if not labeled["rms_omb"] else "")
        labeled["rms_omb"] = True

        plt.barh(y + offsets[3] * bar_width[i], stats.rms_oma[i], height=bar_width[i],
                 color="red", label=labels[3] if not labeled["rms_oma"] else "")
        labeled["rms_oma"] = True

    plt.ylabel("Channel")
    plt.title("Mean & RMS of Obs Residuals vs Channel")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()


def _panel_jo(binned_data: BinnedData, stats: StatisticsData):
    """Panel 3: Jo/p vs Channel."""
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels = ["Jo(b)/p", "Jo(a)/p"]
    labeled = {"job": False, "joa": False}

    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], stats.mean_job[i], height=bar_width[i],
                 color="blue", label=labels[0] if not labeled["job"] else "")
        labeled["job"] = True

        plt.barh(y + offsets[1] * bar_width[i], stats.mean_joa[i], height=bar_width[i],
                 color="red", label=labels[1] if not labeled["joa"] else "")
        labeled["joa"] = True

    plt.ylabel("Channel")
    plt.title("Jo/p vs Channel")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.legend()


def _panel_sigo(binned_data: BinnedData, stats: StatisticsData):
    """Panel 4: Prescribed & Estimated Errors vs Channel."""
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offsets = [-1.5, -0.5, 0.5, 1.5]

    labels = ["sigO", "esigO", "esigB"]
    labeled = {"sigo": False, "esigo": False, "esigb": False}

    for i in range(len(bin_centers)):
        y = bin_centers[i]

        plt.barh(y + offsets[0] * bar_width[i], stats.mean_sigo[i], height=bar_width[i],
                 color="cyan", label=labels[0] if not labeled["sigo"] else "")
        labeled["sigo"] = True

        plt.barh(y + offsets[1] * bar_width[i], stats.mean_esigo[i], height=bar_width[i],
                 color="orange", label=labels[1] if not labeled["esigo"] else "")
        labeled["esigo"] = True

        plt.barh(y + offsets[2] * bar_width[i], stats.mean_esigb[i], height=bar_width[i],
                 color="black", label=labels[2] if not labeled["esigb"] else "")
        labeled["esigb"] = True

    plt.ylabel("Channel")
    plt.title("Prescribed & Estimated Errors vs Channel")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 10))
    plt.legend()






###################################################################################

#utils.py
def timer(base_fn):
    def enhanced_fn():
        start_time = time.perf_counter()
        base_fn()
        end_time = time.perf_counter()
        print(f"Task time: {end_time - start_time} seconds")
    return enhanced_fn    



filename = "python/amsua_metop-b.20260125T150000Z.nc4"
#@timer
def main() -> None:
    
    #IODA file
    reader = IODAReader()
    data = reader.read(filename)

    #masking
    valid_mask = fill_val_mask(data)
    
    #filtering
    valid_data = apply_filter(data, valid_mask)

    #QC masking
    pass_mask = qc_pass_mask(valid_data)
    fail_mask = qc_fail_mask(valid_data)


    pass_data = apply_filter(valid_data, pass_mask)
    fail_data = apply_filter(valid_data, fail_mask)

    #calculate job, joa, esigo, esigb
    pass_data = calc_derived(pass_data)

    #binning
    pass_data_binned = create_channel_bins(pass_data)
    fail_data_binned = create_channel_bins(fail_data)
    #stats
    data_stats = calculate_stats(pass_data_binned)

    #plotting
    panel_plot = plot_stats(pass_data_binned,fail_data_binned, data_stats)
    plt.show()







if __name__ == '__main__':
    main()


