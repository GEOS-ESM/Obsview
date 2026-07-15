#Module for loading IODA data, creating masks, binning levels, and making a 4 panel plot
#TODO: Make this module only have statistics calculation functions
"""Compute (and plot) JEDI-vs-GSI comparison statistics from .nc4 files."""
import numpy as np
import matplotlib.pyplot as plt

from ..plotting.panels import (
    comp_nobs_panel,
    show_resstats_panel,
    show_jo_panel,
    show_sigo_panel,
)
from .binning import get_pressure_mask, accum_pressure_mask

#Function for loading IODA and GSI data, masking missing values and by QC, calculating stats, and making 4 panel plots
#For channel levels
def jediXgsi_channel(nc, varname, levlim, nbins, scaleby, satid, usrqc, comp_bias, common) -> None: #Note the unused parameters 
    #Similar check to varBC in ioda_stats, check if we want to calculate omb or not
    calculate_omb = False   #This is hardcoded for now, needs to become a user selected option later
    obs = nc.groups["ObsValue"].variables[varname][:]
    obs_fill = nc.groups["ObsValue"].variables[varname]._FillValue

    # from JEDI
    #Load relevant data and fill values from JEDI side of the IODA file 
    ombg = nc.groups["ombg"].variables[varname][:]
    ombg_fill = nc.groups["ombg"].variables[varname]._FillValue     #This variable is not used later in the code
    #JEDI H of x
    jhox = nc.groups["hofx0"].variables[varname][:]
    jhox_fill = nc.groups["hofx0"].variables[varname]._FillValue
    #JEDI bias
    jbias = nc.groups["ObsBias0"].variables[varname][:]
    jbias_fill = nc.groups["ObsBias0"].variables[varname]._FillValue
    #JEDI sigo
    jsigo = nc.groups["EffectiveError0"].variables[varname][:]
    jsigo_fill = nc.groups["EffectiveError0"].variables[varname]._FillValue
    #JEDI quality control
    jqc = nc.groups["EffectiveQC0"].variables[varname][:]
    jqc_fill = nc.groups["EffectiveQC0"].variables[varname]._FillValue

    # from GSI:
    #Load relevant data and fill values from the GSI portion of the file
    gomb = nc.groups["GsiHofX"].variables[varname][:]
    gomb_fill = nc.groups["GsiHofX"].variables[varname]._FillValue
    #GSI bias
    gbias = nc.groups["GsiBc"].variables[varname][:]
    gbias_fill = nc.groups["GsiBc"].variables[varname]._FillValue
    #GSI sigo
    gsigo = nc.groups["GsiFinalObsError"].variables[varname][:]
    gsigo_fill = nc.groups["GsiFinalObsError"].variables[varname]._FillValue
    #GSI quality control
    gqc = nc.groups["GsiEffectiveQC"].variables[varname][:]
    gqc_fill = nc.groups["GsiEffectiveQC"].variables[varname]._FillValue

    #Chanel numbers, shared by IODA and GSI
    channel_numbers = nc.groups["MetaData"].variables["sensorChannelNumber"][:]

    #Fill missing values with NaNs or keep data
    obs = np.where(obs == obs_fill, np.nan, obs)
    jsigo = np.where(jsigo == jsigo_fill, np.nan, jsigo)
    gsigo = np.where(gsigo == gsigo_fill, np.nan, gsigo)
    jhox = np.where(jhox == jhox_fill, np.nan, jhox)
    gomb = np.where(gomb == gomb_fill, np.nan, gomb)
    jbias = np.where(jbias == jbias_fill, np.nan, jbias)
    gbias = np.where(gbias == gbias_fill, np.nan, gbias)
    jqc = np.where(jqc == jqc_fill, np.nan, jqc)
    gqc = np.where(gqc == gqc_fill, np.nan, gqc)
    #MISSING: JEDI's omb, used later as 'ombg'
    
    #If parameter 'common' is True,  then gqc and jqc are shared (i.e. only shows the JEDI and GSI output that passes both quality controls)
    if common:
        mask = np.logical_and(gqc == 0, jqc == 0)
        gqc = np.where(mask, gqc, np.nan)
    #NOTE: There isn't the same option for the CLI argument '--bias'

    #Create masked nparray of obs, ombg
    #Put these in the same order as above
    jobs_masked = np.where(jqc == 0, obs, np.nan)
    gobs_masked = np.where(gqc == 0, obs, np.nan)
    ombg_masked = np.where(jqc == 0, ombg, np.nan)
    jhox_masked = np.where(jqc == 0, jhox, np.nan)
    gomb_masked = np.where(gqc == 0, gomb, np.nan)
    jbias_masked = np.where(jqc == 0, jbias, np.nan)
    gbias_masked = np.where(gqc == 0, gbias, np.nan)
    jsigo_masked = np.where(jqc == 0, jsigo, np.nan)
    gsigo_masked = np.where(gqc == 0, gsigo, np.nan)

    #Calculate GSI omb no matter what the calculate option is
    gomb_masked = gobs_masked - gomb_masked - gbias_masked
    #Calculate omb if the option is selected, this option does not exist yet
    if calculate_omb:
        jomb_masked = jobs_masked - jhox_masked - jbias_masked
    else:
        jomb_masked = ombg_masked

    #Calculate stats
    gsum_nobs = np.sum(gqc == 0, axis=0)
    jsum_nobs = np.sum(jqc == 0, axis=0)
    mean_jomb = np.nanmean(jomb_masked, axis=0)
    mean_gomb = np.nanmean(gomb_masked, axis=0)
    mean_jsigo = np.nanmean(jsigo_masked, axis=0)
    mean_gsigo = np.nanmean(gsigo_masked, axis=0)
    rms_jomb = np.sqrt(np.nanmean(jomb_masked**2, axis=0))
    rms_gomb = np.sqrt(np.nanmean(gomb_masked**2, axis=0))

    #Bias stats
    mean_jbias = np.nanmean(jbias_masked, axis=0)
    mean_gbias = np.nanmean(gbias_masked, axis=0)
    rms_jbias = np.sqrt(np.nanmean(jbias_masked**2, axis=0))
    rms_gbias = np.sqrt(np.nanmean(gbias_masked**2, axis=0))

    #Job and joa
    mean_jjob = np.nansum((jomb_masked / jsigo_masked) ** 2, axis=0)
    mean_gjob = np.nansum((gomb_masked / gsigo_masked) ** 2, axis=0)
    mean_jjob = mean_jjob / jsum_nobs
    mean_gjob = mean_gjob / gsum_nobs

    #Create bins for plotting, these are automatic based on number of channels
    bin_centers = np.arange(1, len(channel_numbers) + 1)
    bin_heights = 0.95 * np.ones(len(bin_centers), dtype=int)
    channel_numbers = bin_centers

    #Create 4 panel plot
    plt.figure(figsize=(10, 7))
    plt.subplot(2, 2, 1)
    comp_nobs_panel(varname, jsum_nobs, gsum_nobs, channel_numbers, bin_heights, usrqc, True)
    plt.subplot(2, 2, 2)
    labels = ["JEDI mean", "GSI mean", "JEDI RMS", "GSI RMS"]
    show_resstats_panel(varname, mean_jomb, mean_gomb, rms_jomb, rms_gomb, channel_numbers, bin_heights, labels, True)
    plt.subplot(2, 2, 3)
    labels = ["JEDI Jo(b)/p", "GSI Jo(b)/p"]
    show_jo_panel(varname, mean_jjob, mean_gjob, channel_numbers, bin_heights, labels, True)
    plt.subplot(2, 2, 4)
    #If the compute bias option is selected, the fourth plot will be another resstats panel showing BC stuff
    if comp_bias:
        labels = ["JEDI mean BC", "GSI mean BC", "JEDI RMS BC", "GSI RMS BC"]
        show_resstats_panel(varname, mean_jbias, mean_gbias, rms_jbias, rms_gbias, channel_numbers, bin_heights, labels, True)
    #Otherwise the 4th panel will be the usual sigo and estiamted sigo plots
    else:
        labels = ["JEDI sigO", "GSI sigO"]
        show_sigo_panel(varname, mean_jsigo, mean_gsigo, mean_gsigo, channel_numbers, bin_heights, labels, True)

#Function for loading IODA and GSI data, masking missing values and by QC, binning pressure levels, calculating stats, and making 4 panel plots
#For pressure levels
def jediXgsi_pressure_binned(nc, varname, levlim, nbins, scaleby, satid, usrqc):
    #Load observation values and missing number fill values
    obs = nc.groups["ObsValue"].variables[varname][:]
    obs_fill = nc.groups["ObsValue"].variables[varname]._FillValue

    #Load values into np array and fill missing values with fill value
    # from JEDI:
    ombg = nc.groups["ombg"].variables[varname][:]
    ombg_fill = nc.groups["ombg"].variables[varname]._FillValue
    jsigo = nc.groups["EffectiveError0"].variables[varname][:]
    jsigo_fill = nc.groups["EffectiveError0"].variables[varname]._FillValue
    jqc = nc.groups["EffectiveQC0"].variables[varname][:]
    jqc_fill = nc.groups["EffectiveQC0"].variables[varname]._FillValue

    # from GSI:
    ghox = nc.groups["GsiHofXBc"].variables[varname][:]
    ghox_fill = nc.groups["GsiHofXBc"].variables[varname]._FillValue
    gsigo = nc.groups["GsiFinalObsError"].variables[varname][:]
    gsigo_fill = nc.groups["GsiFinalObsError"].variables[varname]._FillValue
    gqc = nc.groups["GsiEffectiveQC"].variables[varname][:]
    gqc_fill = nc.groups["GsiEffectiveQC"].variables[varname]._FillValue

    #Load pressure levels, shared by JEDI and GSI
    pressure = nc.groups["MetaData"].variables["pressure"][:]
    pressure_fill = nc.groups["MetaData"].variables["pressure"]._FillValue

    #Set satellite identifier variable 'satid' to appropriate value
    if satid > 0:
        sid = nc.groups["MetaData"].variables["satelliteIdentifier"][:]
    else:
        sid = 0

    #Create pressure mask 
    valid_mask = get_pressure_mask(usrqc, sid, pressure, jqc, satid, pressure_fill, jqc_fill)

    #Create bins based on pressure levels that are masked out by the calculated pressure mask
    pressure_valid = pressure[valid_mask] / 100.0
    bins = np.logspace(np.log10(levlim[1]), np.log10(levlim[0]), num=nbins)
    bin_indices = np.digitize(pressure_valid, bins)         #For binning stats
    bin_centers = (bins[:-1] + bins[1:]) / 2                #For plotting
    bin_heights = np.diff(bins)                             #For plotting

    #Create several numpy arrays at once with respective data masked by existing pressure mask
    #pressure mask function also calculates stats
    (jsum_nobs, mean_jomb, rms_jomb, _dum1, _dum2, mean_jjob, _dum3,
     mean_jsigo, _mean_scale) = accum_pressure_mask(
        bins, bin_indices, scaleby, usrqc,
        sid, ombg, ombg, jsigo, pressure, jqc,
        satid, ombg_fill, ombg_fill, jsigo_fill, pressure_fill, jqc_fill,
    )
    #Calculate GSI omb and fill values
    gomb = obs - ghox
    gomb_fill = obs_fill - ghox_fill

    #Create masked stats data
    (gsum_nobs, mean_gomb, rms_gomb, _dum1, _dum2, mean_gjob, _dum3,
     mean_gsigo, _mean_scale) = accum_pressure_mask(
        bins, bin_indices, scaleby, usrqc,
        sid, gomb, ghox, gsigo, pressure, gqc,
        satid, gomb_fill, ghox_fill, gsigo_fill, pressure_fill, gqc_fill,
    )
    #Make 4 panel plot
    plt.figure(figsize=(10, 7))
    plt.subplot(2, 2, 1)
    comp_nobs_panel(varname, jsum_nobs, gsum_nobs, bin_centers, bin_heights, usrqc, False)
    plt.subplot(2, 2, 2)
    labels = ["JEDI mean", "GSI mean", "JEDI RMS", "GSI RMS"]
    show_resstats_panel(varname, mean_jomb, mean_gomb, rms_jomb, rms_gomb, bin_centers, bin_heights, labels, False)
    plt.subplot(2, 2, 3)
    labels = ["JEDI Jo(b)/p", "GSI Jo(b)/p"]
    show_jo_panel(varname, mean_jjob, mean_gjob, bin_centers, bin_heights, labels, False)
    plt.subplot(2, 2, 4)
    labels = ["JEDI sigO", "GSI sigO"]
    show_sigo_panel(varname, mean_jsigo, mean_gsigo, mean_gsigo, bin_centers, bin_heights, labels, False)
