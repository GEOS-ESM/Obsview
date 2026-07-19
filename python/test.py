from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, fields, field, replace
from typing import Optional

filename = "python/x0050.diag_amsua_n19.20230801_00z.ods"


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
    # lev_type: Optional[str] #'pressure' or 'channel'
    #error, also add other relevant variables used in calculations
    #scale
    kt: np.ndarray
    sid: np.ndarray
    #calculated variables
    amb: np.ndarray
    job: np.ndarray
    joa: np.ndarray
    esigo: np.ndarray
    esigb: np.ndarray
    #metadata: Metadata
    all_lev: np.ndarray


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
        job = raw["omb"]**2/raw["sigo"]**2
        joa = raw["oma"]**2/raw["sigo"]**2
        esigo = raw["omb"]*raw["oma"]
        esigb = raw["omb"]*amb

        #Append
        raw["amb"] = amb
        raw["job"] = job
        raw["joa"] = joa
        raw["esigo"] = esigo
        raw["esigb"] = esigb
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
#TODO: Create IODAReader class










###############################
#       Processing Package    #
###############################

#masking.py
def valid_mask(data: ObservationData) -> np.ndarray: 
    missing_val = 1.0e15
    valid_mask = (
        (data.qc == 0)              #Should be changed later to take user input
        #& (data.sid == -999)        #Ditto
        & (data.kt == 40)           #Ditto
        & (data.lev < missing_val)
        & (data.omb < missing_val)
        & (data.oma < missing_val)
        & (data.amb < missing_val)
    )
    return valid_mask





#filtering.py
def apply_filter(data: ObservationData, mask: np.ndarray) -> ObservationData:
    updated_fields = {}
    ignore_fields = {"all_lev"}

    for field in fields(data):
        value = getattr(data, field.name)

        if isinstance(value, np.ndarray) and field.name not in ignore_fields: #If the field is type ndarray and not in the ignore fields dict
            updated_fields[field.name] = value[mask]
        else:
            updated_fields[field.name] = value

    return replace(data, **updated_fields)





#binning.py
@dataclass
class BinnedData:
    data: ObservationData           #Rearranged data by bin
    bin_centers: np.ndarray         #Averaged value between bin levels (used for plotting)
    bin_indices: np.ndarray         
    bin_labels: np.ndarray          #Array of each bin level (unique)
    bin_heights: np.ndarray
    #level_type: str (pressure or channel)

    
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

def plot_stats(binned_data: BinnedData, stats: StatisticsData):
    """
    Build the 4-panel statistics figure from a BinnedData object and a
    StatisticsData object.

    Returns the matplotlib Figure so the caller can later decide to
    plt.show() it or save it with fig.savefig(...).
    """
    fig = plt.figure(figsize=(10, 7))  # hardcoded size for now

    plt.subplot(2, 2, 1)
    _panel_nobs(binned_data, stats)

    plt.subplot(2, 2, 2)
    _panel_resstats(binned_data, stats)

    plt.subplot(2, 2, 3)
    _panel_jo(binned_data, stats)

    plt.subplot(2, 2, 4)
    _panel_sigo(binned_data, stats)

    return fig


def _panel_nobs(binned_data: BinnedData, stats: StatisticsData):
    """Panel 1: Observation count vs Channel."""
    bin_centers = binned_data.bin_centers
    bin_heights = binned_data.bin_heights
    bar_width = bin_heights * 0.4
    offset = -0.5

    labeled = False
    for i in range(len(bin_centers)):
        y = bin_centers[i]
        plt.barh(
            y + offset * bar_width[i],
            stats.nobs[i],
            height=bar_width[i],
            color="green",
            label="used" if not labeled else "",
        )
        labeled = True

    # Channel hardcoded for now (will use binned_data.level_type later).
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

def main() -> None:
    reader = ODSReader()
    data = reader.read(filename)

    print(f"Information about data in file {filename}")
    print(f"Maximum level: {np.max(data.lev)}")
    print(f"Max amb: {np.max(data.amb)}")
    print(f"Number of observations: {np.size(data.obs)}")
    print(f"Length of 'all_lev' array is: {np.size(data.all_lev)}")
    #masking
    my_mask = valid_mask(data)
    print(f"The length of the mask array is: {np.size(my_mask)}")

    #filtering
    filtered_data = apply_filter(data, my_mask)
    print(f"Length of the new observations: {np.size(filtered_data.obs)}")

    #binning
    binned_data = create_channel_bins(filtered_data)

    print(f"Length of of binned data 'obs' array is: {np.size(binned_data.data.obs)}")
    #stats
    stats_data = calculate_stats(binned_data)

    #plotting
    print("Time to plot!")
    my_plot = plot_stats(binned_data, stats_data)
    plt.show()









if __name__ == '__main__':
    main()


