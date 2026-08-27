#Module for creating spatial plots for radiance monitoring
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime

from ..processing.binning import BinnedData
from ..config import KX_TO_DATASRC, KT_TO_DATATYPE

DEFAULT_CHANNEL = 18

def plot_radmon(pass_data: BinnedData, channel: int = DEFAULT_CHANNEL):
    channel_mask = (pass_data.bin_indices == channel)
    values = pass_data.data.omb_no_bias[channel_mask]
    lons = pass_data.data.lon[channel_mask]
    lats = pass_data.data.lat[channel_mask]
    
    
    fig = plt.figure(figsize=(10, 6.5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines(resolution='110m', color='black', linewidth=1)
    ax.set_global()

    sc = ax.scatter(
    lons, lats, 
    c=values, 
    cmap='coolwarm', 
    vmin=-10.0, 
    vmax=10.0, 
    s=3.5,                 # Small marker size to preserve swath visibility
    transform=ccrs.PlateCarree()
    )

    data_src = data_src = KX_TO_DATASRC.get(pass_data.data.kx)
    data_type =  KT_TO_DATATYPE.get(pass_data.data.kt)
    file_type = pass_data.data.file_type
    time = time = pass_data.data.datetime

    title = f"{data_src}\n {time.strftime('%d%b%Y %HZ')}\n All Observations, Channel {channel} 183.310 GHz"

    plt.title(title, fontsize=12, pad=12,multialignment='center')

    cbar = plt.colorbar(sc, orientation='horizontal', pad=0.06, shrink=0.88, aspect=35)

    # Calculate statistics 
    nobs = len(values)
    avg = np.mean(values)
    std = np.std(values)

    cbar.set_label(
    f"Brightness Temperature O-F w/o BC (K)\nnobs={nobs}, avg={avg:.2f}, std={std:.2f}", 
    fontsize=11, 
    labelpad=6
    )
    cbar.ax.tick_params(labelsize=11)
    ...