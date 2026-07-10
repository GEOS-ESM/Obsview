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


def jediXgsi_channel(nc, varname, levlim, nbins, scaleby, satid, usrqc, comp_bias, common):
    calculate_omb = False
    obs = nc.groups["ObsValue"].variables[varname][:]
    obs_fill = nc.groups["ObsValue"].variables[varname]._FillValue

    # from JEDI:
    ombg = nc.groups["ombg"].variables[varname][:]
    ombg_fill = nc.groups["ombg"].variables[varname]._FillValue
    jhox = nc.groups["hofx0"].variables[varname][:]
    jhox_fill = nc.groups["hofx0"].variables[varname]._FillValue
    jbias = nc.groups["ObsBias0"].variables[varname][:]
    jbias_fill = nc.groups["ObsBias0"].variables[varname]._FillValue

    jsigo = nc.groups["EffectiveError0"].variables[varname][:]
    jsigo_fill = nc.groups["EffectiveError0"].variables[varname]._FillValue

    # from GSI:
    gomb = nc.groups["GsiHofX"].variables[varname][:]
    gomb_fill = nc.groups["GsiHofX"].variables[varname]._FillValue
    gbias = nc.groups["GsiBc"].variables[varname][:]
    gbias_fill = nc.groups["GsiBc"].variables[varname]._FillValue

    gsigo = nc.groups["GsiFinalObsError"].variables[varname][:]
    gsigo_fill = nc.groups["GsiFinalObsError"].variables[varname]._FillValue

    gqc = nc.groups["GsiEffectiveQC"].variables[varname][:]
    gqc_fill = nc.groups["GsiEffectiveQC"].variables[varname]._FillValue

    jqc = nc.groups["EffectiveQC0"].variables[varname][:]
    jqc_fill = nc.groups["EffectiveQC0"].variables[varname]._FillValue

    channel_numbers = nc.groups["MetaData"].variables["sensorChannelNumber"][:]

    obs = np.where(obs == obs_fill, np.nan, obs)
    jsigo = np.where(jsigo == jsigo_fill, np.nan, jsigo)
    gsigo = np.where(gsigo == gsigo_fill, np.nan, gsigo)
    jhox = np.where(jhox == jhox_fill, np.nan, jhox)
    gomb = np.where(gomb == gomb_fill, np.nan, gomb)
    jbias = np.where(jbias == jbias_fill, np.nan, jbias)
    gbias = np.where(gbias == gbias_fill, np.nan, gbias)
    jqc = np.where(jqc == jqc_fill, np.nan, jqc)
    gqc = np.where(gqc == gqc_fill, np.nan, gqc)

    if common:
        mask = np.logical_and(gqc == 0, jqc == 0)
        gqc = np.where(mask, gqc, np.nan)

    jobs_masked = np.where(jqc == 0, obs, np.nan)
    gobs_masked = np.where(gqc == 0, obs, np.nan)
    ombg_masked = np.where(jqc == 0, ombg, np.nan)
    jhox_masked = np.where(jqc == 0, jhox, np.nan)
    gomb_masked = np.where(gqc == 0, gomb, np.nan)
    jbias_masked = np.where(jqc == 0, jbias, np.nan)
    gbias_masked = np.where(gqc == 0, gbias, np.nan)
    jsigo_masked = np.where(jqc == 0, jsigo, np.nan)
    gsigo_masked = np.where(gqc == 0, gsigo, np.nan)

    gomb_masked = gobs_masked - gomb_masked - gbias_masked
    if calculate_omb:
        jomb_masked = jobs_masked - jhox_masked - jbias_masked
    else:
        jomb_masked = ombg_masked

    gsum_nobs = np.sum(gqc == 0, axis=0)
    jsum_nobs = np.sum(jqc == 0, axis=0)
    mean_jomb = np.nanmean(jomb_masked, axis=0)
    mean_gomb = np.nanmean(gomb_masked, axis=0)
    mean_jsigo = np.nanmean(jsigo_masked, axis=0)
    mean_gsigo = np.nanmean(gsigo_masked, axis=0)
    rms_jomb = np.sqrt(np.nanmean(jomb_masked**2, axis=0))
    rms_gomb = np.sqrt(np.nanmean(gomb_masked**2, axis=0))

    mean_jbias = np.nanmean(jbias_masked, axis=0)
    mean_gbias = np.nanmean(gbias_masked, axis=0)
    rms_jbias = np.sqrt(np.nanmean(jbias_masked**2, axis=0))
    rms_gbias = np.sqrt(np.nanmean(gbias_masked**2, axis=0))

    mean_jjob = np.nansum((jomb_masked / jsigo_masked) ** 2, axis=0)
    mean_gjob = np.nansum((gomb_masked / gsigo_masked) ** 2, axis=0)
    mean_jjob = mean_jjob / jsum_nobs
    mean_gjob = mean_gjob / gsum_nobs

    bin_centers = np.arange(1, len(channel_numbers) + 1)
    bin_heights = 0.95 * np.ones(len(bin_centers), dtype=int)
    channel_numbers = bin_centers

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
    if comp_bias:
        labels = ["JEDI mean BC", "GSI mean BC", "JEDI RMS BC", "GSI RMS BC"]
        show_resstats_panel(varname, mean_jbias, mean_gbias, rms_jbias, rms_gbias, channel_numbers, bin_heights, labels, True)
    else:
        labels = ["JEDI sigO", "GSI sigO"]
        show_sigo_panel(varname, mean_jsigo, mean_gsigo, mean_gsigo, channel_numbers, bin_heights, labels, True)


def jediXgsi_pressure_binned(nc, varname, levlim, nbins, scaleby, satid, usrqc):
    obs = nc.groups["ObsValue"].variables[varname][:]
    obs_fill = nc.groups["ObsValue"].variables[varname]._FillValue

    # from JEDI:
    ombg = nc.groups["ombg"].variables[varname][:]
    ombg_fill = nc.groups["ombg"].variables[varname]._FillValue
    jsigo = nc.groups["EffectiveError0"].variables[varname][:]
    jsigo_fill = nc.groups["EffectiveError0"].variables[varname]._FillValue

    # from GSI:
    ghox = nc.groups["GsiHofXBc"].variables[varname][:]
    ghox_fill = nc.groups["GsiHofXBc"].variables[varname]._FillValue
    gsigo = nc.groups["GsiFinalObsError"].variables[varname][:]
    gsigo_fill = nc.groups["GsiFinalObsError"].variables[varname]._FillValue

    gqc = nc.groups["GsiEffectiveQC"].variables[varname][:]
    gqc_fill = nc.groups["GsiEffectiveQC"].variables[varname]._FillValue

    jqc = nc.groups["EffectiveQC0"].variables[varname][:]
    jqc_fill = nc.groups["EffectiveQC0"].variables[varname]._FillValue

    pressure = nc.groups["MetaData"].variables["pressure"][:]
    pressure_fill = nc.groups["MetaData"].variables["pressure"]._FillValue

    if satid > 0:
        sid = nc.groups["MetaData"].variables["satelliteIdentifier"][:]
    else:
        sid = 0

    valid_mask = get_pressure_mask(usrqc, sid, pressure, jqc, satid, pressure_fill, jqc_fill)

    pressure_valid = pressure[valid_mask] / 100.0
    bins = np.logspace(np.log10(levlim[1]), np.log10(levlim[0]), num=nbins)
    bin_indices = np.digitize(pressure_valid, bins)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_heights = np.diff(bins)

    (jsum_nobs, mean_jomb, rms_jomb, _dum1, _dum2, mean_jjob, _dum3,
     mean_jsigo, _mean_scale) = accum_pressure_mask(
        bins, bin_indices, scaleby, usrqc,
        sid, ombg, ombg, jsigo, pressure, jqc,
        satid, ombg_fill, ombg_fill, jsigo_fill, pressure_fill, jqc_fill,
    )

    gomb = obs - ghox
    gomb_fill = obs_fill - ghox_fill

    (gsum_nobs, mean_gomb, rms_gomb, _dum1, _dum2, mean_gjob, _dum3,
     mean_gsigo, _mean_scale) = accum_pressure_mask(
        bins, bin_indices, scaleby, usrqc,
        sid, gomb, ghox, gsigo, pressure, gqc,
        satid, gomb_fill, ghox_fill, gsigo_fill, pressure_fill, gqc_fill,
    )

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
