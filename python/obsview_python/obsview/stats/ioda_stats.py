#Module for creating binned statistics and plots for IODA files
"""Compute (and plot) statistics from IODA/JEDI .nc4 files."""
import numpy as np
import matplotlib.pyplot as plt

from ..plotting.panels import (
    show_nobs_panel,
    show_resstats_panel,
    show_jo_panel,
    show_sigo_panel,
)

#Function for loading IODA data, masking missing values, creating bins, calculating stats, creating 4 panel plot
def ioda_pressure_binned(nc, varname, levlim, nbins, scaleby, satid, usrqc):
    #Load data
    ombg = nc.groups["ombg"].variables[varname][:]
    oman = nc.groups["oman"].variables[varname][:]
    sigo = nc.groups["EffectiveError0"].variables[varname][:]
    #If scale factor specified, add scale data
    if scaleby != "null":
        scale = nc.groups[scaleby].variables[varname][:]
    qc = nc.groups["EffectiveQC0"].variables[varname][:]
    pressure = nc.groups["MetaData"].variables["pressure"][:]
    if satid > 0:
        sid = nc.groups["MetaData"].variables["satelliteIdentifier"][:]
    amb = ombg - oman #Calculate analysis minus observation

    #Add fill values
    ombg_fill = nc.groups["ombg"].variables[varname]._FillValue
    oman_fill = nc.groups["oman"].variables[varname]._FillValue
    sigo_fill = nc.groups["EffectiveError0"].variables[varname]._FillValue
    #If scale factor included, add scale fill value
    if scaleby != "null":
        scale_fill = nc.groups[scaleby].variables[varname]._FillValue
    pressure_fill = nc.groups["MetaData"].variables["pressure"]._FillValue
    qc_fill = nc.groups["EffectiveQC0"].variables[varname]._FillValue

    #Create masks for missing value data
    if satid > 0:   #If satid specified, then the mask includes satid booleans
        valid_mask = (
            (sid == satid)
            & (qc == usrqc)
            & (ombg != ombg_fill)
            & (oman != oman_fill)
            & (sigo != sigo_fill)
            & (pressure != pressure_fill)
            & (qc != qc_fill)
        )
    else: #Otherwise the mask will check if relevant variables are ALL true
        valid_mask = (
            (qc == usrqc)
            & (ombg != ombg_fill)
            & (oman != oman_fill)
            & (sigo != sigo_fill)
            & (pressure != pressure_fill)
            & (qc != qc_fill)
        )

    #Add mask to the data
    nobs_valid = qc[valid_mask]
    ombg_valid = ombg[valid_mask]
    oman_valid = oman[valid_mask]
    sigo_valid = sigo[valid_mask]
    amb_valid = amb[valid_mask]
    if scaleby != "null":
        scale_valid = scale[valid_mask]

    #Create pressure bins
    pressure_valid = pressure[valid_mask] / 100.0
    bins = np.logspace(np.log10(levlim[1]), np.log10(levlim[0]), num=nbins)
    bin_indices = np.digitize(pressure_valid, bins)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_heights = np.diff(bins)

    #Initialize stats arrays
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

    #Loop over each bin
    for i in range(1, len(bins)):
        bin_mask = bin_indices == i
        if np.any(bin_mask): #If the bin mask is true...
            #Bin relevant data
            ombg_bin = ombg_valid[bin_mask]
            oman_bin = oman_valid[bin_mask]
            amb_bin = amb_valid[bin_mask]
            sigo_bin = sigo_valid[bin_mask]

            #calculate stats and append them to their respective lists
            sum_nobs.append(int(np.sum(bin_mask)))          #Total number of observations
            mean_ombg.append(np.mean(ombg_bin))             #Mean omb
            rms_ombg.append(np.sqrt(np.mean(ombg_bin**2)))  #RMS omb
            mean_oman.append(np.mean(oman_bin))             #Mean oma
            rms_oman.append(np.sqrt(np.mean(oman_bin**2)))  #RMS oma

            #Calculate actual mean job and joa
            mean_job.append(np.sum((ombg_bin / sigo_bin) ** 2))
            mean_joa.append(np.sum((oman_bin / sigo_bin) ** 2))

            #Calculate mean sigo and estimated mean sig o
            mean_sigo.append(np.mean(sigo_bin))
            mean_esigo.append(np.sqrt(np.abs(np.mean(ombg_bin * oman_bin))))
            mean_esigb.append(np.sqrt(np.abs(np.mean(ombg_bin * amb_bin))))

            if scaleby != "null":
                scale_bin = scale_valid[bin_mask]
                mean_scale.append(np.mean(scale_bin))
        else: #If bin mask is false, append NaN to stats list row
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

    #Scale stats by pressure level if scaling is chosen as an option 
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


#Function for loading IODA data, masking missing values, creating bins, calculating stats, creating 4 panel plot
def ioda_channel(nc, varname, levlim, nbins, scaleby, satid, usrqc): #Note the parameters that are not used
    #Variational bias correction: changes where and how omb(and others) is sourced
    varBC = False #Hardcoded that we trust omb values directly from IODA file
    if varBC:
        obs = nc.groups["ObsValue"].variables[varname][:]
        obs_fill = nc.groups["ObsValue"].variables[varname]._FillValue
        ombg = nc.groups["hofx0"].variables[varname][:]
        ombg_fill = nc.groups["hofx0"].variables[varname]._FillValue
        bias = nc.groups["ObsBias0"].variables[varname][:]
        bias_fill = nc.groups["ObsBias0"].variables[varname]._FillValue
        ombg = obs - ombg + bias    #Manually calculate omb
        ombg_fill = obs_fill - ombg_fill + bias_fill
    else:
        ombg = nc.groups["ombg"].variables[varname][:]
        ombg_fill = nc.groups["ombg"].variables[varname]._FillValue
    #Load other necessary data and their fill values
    oman = nc.groups["oman"].variables[varname][:]
    oman_fill = nc.groups["oman"].variables[varname]._FillValue
    sigo = nc.groups["EffectiveError0"].variables[varname][:]
    sigo_fill = nc.groups["EffectiveError0"].variables[varname]._FillValue

    #Calculate analysis minus background and its fill values
    amb = ombg - oman
    amb_fill = ombg_fill - oman_fill

    #Quality control and fill values
    qc = nc.groups["EffectiveQC0"].variables[varname][:]
    qc_fill = nc.groups["EffectiveQC0"].variables[varname]._FillValue

    #Where to get channel numbers(everything except aero uses the Channel attritbute)
    if varname == "aerosolOpticalDepth":
        channel_numbers = nc.groups["MetaData"].variables["obs_wavelength"][:]
    else:
        channel_numbers = nc.variables["Channel"][:]

    #Fill missing values with NaNs or data
    ombg = np.where(ombg == ombg_fill, np.nan, ombg)
    oman = np.where(oman == oman_fill, np.nan, oman)
    sigo = np.where(sigo == sigo_fill, np.nan, sigo)
    amb = np.where(amb == amb_fill, np.nan, amb)
    qc = np.where(qc == qc_fill, np.nan, qc)

    #Create quality control masks and keep only data where qc=0, everything else is NaNs
    ombg_masked = np.where(qc == 0, ombg, np.nan)
    oman_masked = np.where(qc == 0, oman, np.nan)
    amb_masked = np.where(qc == 0, amb, np.nan)
    sigo_masked = np.where(qc == 0, sigo, np.nan)

    #Calculate statistics ignoring NaN values
    sum_nobs = np.sum(qc == 0, axis=0)
    mean_ombg = np.nanmean(ombg_masked, axis=0)
    mean_oman = np.nanmean(oman_masked, axis=0)
    mean_sigo = np.nanmean(sigo_masked, axis=0)
    rms_ombg = np.sqrt(np.nanmean(ombg_masked**2, axis=0))
    rms_oman = np.sqrt(np.nanmean(oman_masked**2, axis=0))

    #Calculate estimated sig o
    mean_esigo = np.sqrt(np.abs(np.nanmean(ombg_masked * oman_masked, axis=0)))
    mean_esigb = np.sqrt(np.abs(np.nanmean(ombg_masked * amb_masked, axis=0)))

    #Calculate mean job and joa
    mean_job = np.nansum((ombg_masked / sigo_masked) ** 2, axis=0)
    mean_joa = np.nansum((oman_masked / sigo_masked) ** 2, axis=0)
    mean_job = mean_job / sum_nobs
    mean_joa = mean_joa / sum_nobs

    #Create bins for plotting (bins come directly from channel data so not computation is required)
    bin_centers = np.arange(1, len(channel_numbers) + 1)
    bin_heights = 0.95 * np.ones(len(bin_centers), dtype=int)
    channel_numbers = bin_centers

    #Create 4 panel plots
    plt.figure(figsize=(10, 7))
    plt.subplot(2, 2, 1)
    show_nobs_panel(varname, sum_nobs, sum_nobs, channel_numbers, bin_heights, usrqc, False, True)
    plt.subplot(2, 2, 2)
    labels = ["Mean o-b", "Mean o-a", "RMS o-b", "RMS o-a"]
    show_resstats_panel(varname, mean_ombg, mean_oman, rms_ombg, rms_oman, channel_numbers, bin_heights, labels, True)
    plt.subplot(2, 2, 3)
    labels = ["Jo(b)/p", "Jo(a)/p"]
    show_jo_panel(varname, mean_job, mean_joa, bin_centers, bin_heights, labels, True)
    plt.subplot(2, 2, 4)
    labels = ["sigO", "esigO", "esigB"]
    show_sigo_panel(varname, mean_sigo, mean_esigo, mean_esigb, channel_numbers, bin_heights, labels, True)
