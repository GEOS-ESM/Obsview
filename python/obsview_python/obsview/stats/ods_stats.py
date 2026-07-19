#Module for computing statistics, plotting them, and reading ODS files???
#TODO: Rewrite ods_pressure_binned() and ods_channel() with separate functions for each purpose
"""Compute (and plot) statistics from .ods diagnostic files.""" #Change this
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import xarray as xr
import pandas as pd

from ..plotting.panels import (
    show_nobs_panel,
    show_resstats_panel,
    show_jo_panel,
    show_sigo_panel,
)

#Function to read data binned by pressure levels (not radiance)
def ods_pressure_binned(filename, varname, levlim, nbins, scaleby, satid, kt, usrqc):
    #Open ODS file using netCDF4 module
    with Dataset(filename, "r") as nc:
        #Leaves missing values in array
        nc.set_auto_mask(False)
        #Assigns data to respective variables
        ktp = nc.variables["kt"][:].astype(np.int16)
        sid = nc.variables["kx"][:].astype(np.int32)
        qc = nc.variables["qcexcl"][:]
        pressure = nc.variables["lev"][:]
        ombg = nc.variables["omf"][:]
        oman = nc.variables["oma"][:]
        obs = nc.variables["obs"][:]
        sigo = nc.variables["xvec"][:]
        #Flattens multi-dimensional arrays into one dimension(row-major)
        pressure = pressure.flatten()
        ktp = ktp.flatten()
        sid = sid.flatten()
        qc = qc.flatten()
        ombg = ombg.flatten()
        oman = oman.flatten()
        scale = obs.flatten()
        sigo = sigo.flatten()
        amb = ombg - oman
       #Set scale size if user wanted scaled levels 
        if scaleby == "hofx0" or scaleby == "bkg":
            scale = scale - ombg
    #Create mask to check for missing values and leave clean arrays 
    missing_val = 1.0e15
    valid_mask = (
        (qc == usrqc)
        & (sid == satid)
        & (ktp == kt)
        & (pressure < missing_val)
        & (ombg < missing_val)
        & (oman < missing_val)
        & (amb < missing_val)
    )

    nobs_valid = qc[valid_mask]
    ombg_valid = ombg[valid_mask]
    oman_valid = oman[valid_mask]
    amb_valid = amb[valid_mask]
    sigo_valid = sigo[valid_mask]
    pressure_valid = pressure[valid_mask]
    
    if scaleby != "null":
        scale_valid = scale[valid_mask]
    

    #Creates bins in log space(how pressure is leveled)
    bins = np.logspace(np.log10(levlim[1]), np.log10(levlim[0]), num=nbins)
    bin_indices = np.digitize(pressure_valid, bins)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_heights = np.diff(bins)
    #Initializes arrays for ODS stats objects
    #TODO: make an ODS stats class to handle more cleanly loading this information
    # (ideally using numpy arrays instead of lists)
    sum_nobs = []
    mean_ombg = []
    rms_ombg = []
    mean_oman = []
    rms_oman = []
    mean_job = []
    mean_joa = []
    mean_sigo = []
    mean_esigo = []
    mean_esigb = []
    if scaleby != "null":
        mean_scale = []
    #Loop over bin indices
    for i in range(1, len(bins)):
        bin_mask = bin_indices == i
        #If bins contain data (not all missing values)
        if np.any(bin_mask):
            ombg_bin = ombg_valid[bin_mask]
            oman_bin = oman_valid[bin_mask]
            amb_bin = amb_valid[bin_mask]
            sigo_bin = sigo_valid[bin_mask]

            sum_nobs.append(int(np.sum(bin_mask)))
            mean_ombg.append(np.mean(ombg_bin))
            rms_ombg.append(np.sqrt(np.mean(ombg_bin**2)))
            mean_oman.append(np.mean(oman_bin))
            rms_oman.append(np.sqrt(np.mean(oman_bin**2)))

            mean_job.append(np.sum((ombg_bin / sigo_bin) ** 2))
            mean_joa.append(np.sum((oman_bin / sigo_bin) ** 2))

            mean_sigo.append(np.mean(sigo_bin))
            mean_esigo.append(np.sqrt(np.abs(np.mean(ombg_bin * oman_bin))))
            mean_esigb.append(np.sqrt(np.abs(np.mean(ombg_bin * amb_bin))))

            if scaleby != "null":
                scale_bin = scale_valid[bin_mask]
                mean_scale.append(np.mean(scale_bin))
        #If bin is empty, then fill stats lists with NaNs
        else:
            sum_nobs.append(np.nan)
            mean_ombg.append(np.nan)
            rms_ombg.append(np.nan)
            mean_oman.append(np.nan)
            rms_oman.append(np.nan)
            mean_job.append(np.nan)
            mean_joa.append(np.nan)
            mean_sigo.append(np.nan)
            mean_esigo.append(np.nan)
            mean_esigb.append(np.nan)
            if scaleby != "null":
                mean_scale.append(np.nan)

    #Turn stats lists into numpy arrays
    sum_nobs = np.array(sum_nobs)
    mean_ombg = np.array(mean_ombg)
    rms_ombg = np.array(rms_ombg)
    mean_oman = np.array(mean_oman)
    rms_oman = np.array(rms_oman)
    mean_job = np.array(mean_job)
    mean_joa = np.array(mean_joa)
    mean_sigo = np.array(mean_sigo)
    mean_esigo = np.array(mean_esigo)
    mean_esigb = np.array(mean_esigb)

    #If there is a scaling factor, scale respective stats in each bin
    if scaleby != "null":
        mean_scale = np.array(mean_scale)
        mean_ombg = mean_ombg / mean_scale
        rms_ombg = rms_ombg / mean_scale
        mean_oman = mean_oman / mean_scale
        rms_oman = rms_oman / mean_scale
        mean_sigo = mean_sigo / mean_scale
        mean_esigo = mean_esigo / mean_scale
        mean_esigb = mean_esigb / mean_scale

    #Calculate mean job and joa
    mean_job = mean_job / sum_nobs
    mean_joa = mean_joa / sum_nobs

    #Create 4 panel plot
    plt.figure(figsize=(10, 7))
    plt.subplot(2, 2, 1)
    show_nobs_panel(varname, sum_nobs, sum_nobs, bin_centers, bin_heights, usrqc, False, False)
    plt.subplot(2, 2, 2)
    labels = ["Mean o-b", "Mean o-a", "RMS o-b", "RMS o-a"]
    show_resstats_panel(varname, mean_ombg, mean_oman, rms_ombg, rms_oman, bin_centers, bin_heights, labels, False)
    plt.subplot(2, 2, 3)
    labels = ["Jo(b)/p", "Jo(a)/p"]
    show_jo_panel(varname, mean_job, mean_joa, bin_centers, bin_heights, labels, False)
    plt.subplot(2, 2, 4)
    labels = ["sigO", "esigO", "esigB"]
    show_sigo_panel(varname, mean_sigo, mean_esigo, mean_esigb, bin_centers, bin_heights, labels, False)

#Function to read ODS files, load data, calculate stats, and plot 4 panel figure for data by channel level instead of pressure
def ods_channel(filename, varname, levlim, nbins, satid, kt, usrqc):
    #Open file using xarray
    ds = xr.open_dataset(filename)

    #Load data
    lev = ds["lev"]
    omf = ds["omf"]
    oma = ds["oma"]
    sigo = ds["xvec"]
    qcexcl = ds["qcexcl"]

    #Create single dimensional array
    lev_flat = lev.values.flatten()
    omf_flat = omf.values.flatten()
    oma_flat = oma.values.flatten()
    sigo_flat = sigo.values.flatten()
    qcexcl_flat = qcexcl.values.flatten()

    #Calculate analysis minus background
    amb_flat = omf_flat - oma_flat
    #Estimated sigo and sigb
    esigo_flat = omf_flat * oma_flat
    esigb_flat = omf_flat * amb_flat
    #Calculate JOB and JOA
    job = omf_flat * omf_flat / (sigo_flat * sigo_flat)
    joa = oma_flat * oma_flat / (sigo_flat * sigo_flat)

    #Store flattened stats and respective levels in pandas Dataframe class
    df_all = pd.DataFrame(
        {
            "lev": lev_flat,
            "omf": omf_flat,
            "oma": oma_flat,
            "sigo": sigo_flat,
            "esigo": esigo_flat,
            "esigb": esigb_flat,
            "job": job,
            "joa": joa,
            "qcexcl": qcexcl_flat,
        }
    ).dropna(subset=["lev", "qcexcl"]) #and remove missing values from levels and quality control data
    #Identify all unique atmospheric levels
    all_levs = np.unique(df_all["lev"])

    #Isolate quality controlled data
    df_valid = df_all[df_all["qcexcl"] == 0].dropna(subset=["omf", "oma"])

    #Group dataset by atmospheric level and perform per level statistical aggregation
    grouped = df_valid.groupby("lev").agg(
        {
            "omf": ["mean", lambda x: np.sqrt(np.mean(x**2))],
            "oma": ["mean", lambda x: np.sqrt(np.mean(x**2))],
            "sigo": ["mean", lambda x: np.sqrt(np.mean(x**2))],
            "esigo": [lambda x: np.sqrt(np.abs(np.mean(x)))],
            "esigb": [lambda x: np.sqrt(np.abs(np.mean(x)))],
            "job": ["sum"],
            "joa": ["sum"],
        }
    )
    #Flatten and rename multi-index columns
    grouped.columns = [
        "omf_mean", "omf_rms", "oma_mean", "oma_rms", "sigo_mean", "sigo_rms",
        "esigo_mean", "esigb_mean", "sum_job", "sum_joa",
    ]
    #Convert grouping key back into column
    grouped = grouped.reset_index()

    #Construct complete profile results, including empty levels into result DataFrame object
    full_result = pd.DataFrame({"lev": all_levs})
    result = pd.merge(full_result, grouped, on="lev", how="left").sort_values("lev")

    #Adds additional column that identifies if row passes or fails quality control 
    df_all["qcexcl_valid"] = df_all["qcexcl"] == 0
    #Group data by level and vailidity status
    count_data = df_all.groupby(["lev", "qcexcl_valid"]).size().unstack(fill_value=0)

    #Substitutes header boolean values with descriptive names
    count_data = count_data.rename(columns={True: "valid_count", False: "excluded_count"}).reset_index()

    #Construct final comprehensive volumetric count table
    count_full = pd.merge(full_result, count_data, on="lev", how="left").fillna(0)
    count_full = count_full.sort_values("lev")

    #Create bin centers and heights
    channel_numbers = result["lev"]
    bins = np.arange(len(channel_numbers) + 1)
    bin_centers = np.arange(1, len(channel_numbers) + 1)
    bin_heights = 0.95 * np.ones(len(bin_centers), dtype=int)
    channel_numbers = bin_centers

    #Calculate stats
    sum_nobs = count_full["valid_count"]        #Total number of observations
    sum_job = result["sum_job"] / sum_nobs      #Mean job
    sum_joa = result["sum_joa"] / sum_nobs      #Mean joa


    #Create 4 panel plots
    #TODO: Make this a separate function
    plt.figure(figsize=(10, 7))
    plt.subplot(2, 2, 1)
    show_nobs_panel(varname, sum_nobs, count_full["excluded_count"], bin_centers, bin_heights, usrqc, False, True)
    plt.subplot(2, 2, 2)
    labels = ["Mean o-b", "Mean o-a", "RMS o-b", "RMS o-a"]
    show_resstats_panel(
        varname, result["omf_mean"], result["oma_mean"], result["omf_rms"], result["oma_rms"],
        bin_centers, bin_heights, labels, True,
    )
    plt.subplot(2, 2, 3)
    labels = ["Jo(b)/p", "Jo(a)/p"]
    show_jo_panel(varname, sum_job, sum_joa, bin_centers, bin_heights, labels, True)
    plt.subplot(2, 2, 4)
    labels = ["sigO", "esigO", "esigB"]
    show_sigo_panel(
        varname, result["sigo_mean"], result["esigo_mean"], result["esigb_mean"],
        bin_centers, bin_heights, labels, True,
    )
