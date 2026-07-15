#Module for functions that create bins and make masks for both pressure and channel levels

"""Shared bin-construction and masking/accumulation helpers.

These were originally duplicated (with small variations) inside each of
the pressure-binned stats functions. ``get_pressure_mask`` /
``accum_pressure_mask`` are used by the JEDI-vs-GSI comparison path;
the other pressure-binned functions inline an equivalent mask+accumulate
loop today (see stats/ods_stats.py and stats/ioda_stats.py) and are good
candidates to consolidate onto these helpers in a follow-up pass.
"""
import numpy as np

#Make bin indices, bin centers, and bin heights from pressure levels
def make_pressure_bins(levlim, nbins):
    """Build log-spaced pressure bin edges/centers/heights.

    Parameters
    ----------
    levlim : (bottom, top) in hPa
    nbins : number of bin edges

    Returns
    -------
    bins, bin_centers, bin_heights
    """
    bins = np.logspace(np.log10(levlim[1]), np.log10(levlim[0]), num=nbins)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_heights = np.diff(bins)
    return bins, bin_centers, bin_heights

#Make bin indices, bin centers, and bin heights from channel levels
def make_channel_bins(channel_numbers):
    """Build unit-spaced 'bins' for channel/level-indexed (non-pressure) data."""
    bin_centers = np.arange(1, len(channel_numbers) + 1)
    bin_heights = 0.95 * np.ones(len(bin_centers), dtype=int)
    return bin_centers, bin_centers, bin_heights


#Make a pressure mask
def get_pressure_mask(usrqc, sid, pressure, qc, satid, pressure_fill, qc_fill):
    """Boolean validity mask for pressure-binned data (no ombg/oman/sigo check)."""
    if satid > 0:
        return (
            (sid == satid)
            & (qc == usrqc)
            & (pressure != pressure_fill)
            & (qc != qc_fill)
        )
    return (qc == usrqc) & (pressure != pressure_fill) & (qc != qc_fill)

#Function that masks respective variables by their missing values AND calculated necessary stats
def accum_pressure_mask(            #Incredible amount of parameters, need to reduce this with specific xGSI class
    bins,
    bin_indices,
    scaleby,
    usrqc,
    sid,
    ombg,
    oman,
    sigo,
    pressure,
    qc,
    satid,
    ombg_fill,
    oman_fill,
    sigo_fill,
    pressure_fill,
    qc_fill,
):
    """Mask by QC/fill-value validity, then bin and accumulate o-b/o-a stats.

    Returns arrays (one entry per bin, in the order the bins were given):
    sum_nobs, mean_ombg, rms_ombg, mean_oman, rms_oman, mean_job, mean_joa,
    mean_sigo, mean_scale (mean_scale is an empty array when scaleby=='null').
    """
    #If the satellite ID set to something, then include it in the masking
    if satid > 0:
        valid_mask = (
            (sid == satid)
            & (qc == usrqc)
            & (ombg != ombg_fill)
            & (oman != oman_fill)
            & (sigo != sigo_fill)
            & (pressure != pressure_fill)
            & (qc != qc_fill)
        )
    else:       #Otherwise create masks with the regular variables
        valid_mask = (
            (qc == usrqc)
            & (ombg != ombg_fill)
            & (oman != oman_fill)
            & (sigo != sigo_fill)
            & (pressure != pressure_fill)
            & (qc != qc_fill)
        )

    #Slice the respective data
    nobs_valid = qc[valid_mask]
    ombg_valid = ombg[valid_mask]
    oman_valid = oman[valid_mask]
    sigo_valid = sigo[valid_mask]
    scale = scale - ombg #From ods_stats
    if scaleby != "null":
        scale_valid = scale[valid_mask]  # noqa: F821 - matches original behavior

    #Initialize calculated stats lists
    sum_nobs, mean_ombg, rms_ombg = [], [], []
    mean_oman, rms_oman = [], []
    mean_job, mean_joa = [], []
    mean_sigo = []
    mean_scale = []

    #Loop through each bin
    for i in range(1, len(bins)):
        bin_mask = bin_indices == i
        if np.any(bin_mask):    
            #First mask only the variables of a specific bin
            nobs_bin = nobs_valid[bin_mask]
            ombg_bin = ombg_valid[bin_mask]
            oman_bin = oman_valid[bin_mask]
            sigo_bin = sigo_valid[bin_mask]
            
            #Append calculated stats to their respective lists
            sum_nobs.append(len(nobs_bin))
            mean_ombg.append(np.mean(ombg_bin))
            rms_ombg.append(np.sqrt(np.mean(ombg_bin**2)))
            mean_oman.append(np.mean(oman_bin))
            rms_oman.append(np.sqrt(np.mean(oman_bin**2)))
            mean_job.append(np.sum((ombg_bin / sigo_bin) ** 2))
            mean_joa.append(np.sum((oman_bin / sigo_bin) ** 2))
            mean_sigo.append(np.mean(sigo_bin))

            #If there is a scale factor, then add 
            if scaleby != "null":
                scale_bin = scale_valid[bin_mask]
                mean_scale.append(np.mean(scale_bin))
        else:   #If the bin index is a missing value, then put NaNs for all stats
            sum_nobs.append(np.nan)
            mean_ombg.append(np.nan)
            rms_ombg.append(np.nan)
            mean_oman.append(np.nan)
            rms_oman.append(np.nan)
            mean_job.append(np.nan)
            mean_joa.append(np.nan)
            mean_sigo.append(np.nan)
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

    #Adjust stats by scale factor if user specifies they want things scaled
    if scaleby != "null":
        mean_scale = np.array(mean_scale)
        mean_ombg = mean_ombg / mean_scale
        rms_ombg = rms_ombg / mean_scale
        mean_oman = mean_oman / mean_scale
        rms_oman = rms_oman / mean_scale
        mean_sigo = mean_sigo / mean_scale

    #Calculate mean job and joa
    mean_job = mean_job / sum_nobs
    mean_joa = mean_joa / sum_nobs

    #Return numpy arrays of the calcualted stats and scale factor
    return (            #Incredible amount of returned arrays, make this cleaner by returning an object of CalcdStats class
        sum_nobs,
        mean_ombg,
        rms_ombg,
        mean_oman,
        rms_oman,
        mean_job,
        mean_joa,
        mean_sigo,
        mean_scale,
    )
