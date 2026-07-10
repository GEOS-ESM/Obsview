"""Shared bin-construction and masking/accumulation helpers.

These were originally duplicated (with small variations) inside each of
the pressure-binned stats functions. ``get_pressure_mask`` /
``accum_pressure_mask`` are used by the JEDI-vs-GSI comparison path;
the other pressure-binned functions inline an equivalent mask+accumulate
loop today (see stats/ods_stats.py and stats/ioda_stats.py) and are good
candidates to consolidate onto these helpers in a follow-up pass.
"""
import numpy as np


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


def make_channel_bins(channel_numbers):
    """Build unit-spaced 'bins' for channel/level-indexed (non-pressure) data."""
    bin_centers = np.arange(1, len(channel_numbers) + 1)
    bin_heights = 0.95 * np.ones(len(bin_centers), dtype=int)
    return bin_centers, bin_centers, bin_heights


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


def accum_pressure_mask(
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
    else:
        valid_mask = (
            (qc == usrqc)
            & (ombg != ombg_fill)
            & (oman != oman_fill)
            & (sigo != sigo_fill)
            & (pressure != pressure_fill)
            & (qc != qc_fill)
        )

    nobs_valid = qc[valid_mask]
    ombg_valid = ombg[valid_mask]
    oman_valid = oman[valid_mask]
    sigo_valid = sigo[valid_mask]
    if scaleby != "null":
        scale_valid = scale[valid_mask]  # noqa: F821 - matches original behavior

    sum_nobs, mean_ombg, rms_ombg = [], [], []
    mean_oman, rms_oman = [], []
    mean_job, mean_joa = [], []
    mean_sigo = []
    mean_scale = []

    for i in range(1, len(bins)):
        bin_mask = bin_indices == i
        if np.any(bin_mask):
            nobs_bin = nobs_valid[bin_mask]
            ombg_bin = ombg_valid[bin_mask]
            oman_bin = oman_valid[bin_mask]
            sigo_bin = sigo_valid[bin_mask]

            sum_nobs.append(len(nobs_bin))
            mean_ombg.append(np.mean(ombg_bin))
            rms_ombg.append(np.sqrt(np.mean(ombg_bin**2)))
            mean_oman.append(np.mean(oman_bin))
            rms_oman.append(np.sqrt(np.mean(oman_bin**2)))

            mean_job.append(np.sum((ombg_bin / sigo_bin) ** 2))
            mean_joa.append(np.sum((oman_bin / sigo_bin) ** 2))

            mean_sigo.append(np.mean(sigo_bin))

            if scaleby != "null":
                scale_bin = scale_valid[bin_mask]
                mean_scale.append(np.mean(scale_bin))
        else:
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

    sum_nobs = np.array(sum_nobs)
    mean_ombg = np.array(mean_ombg)
    rms_ombg = np.array(rms_ombg)
    mean_oman = np.array(mean_oman)
    rms_oman = np.array(rms_oman)
    mean_job = np.array(mean_job)
    mean_joa = np.array(mean_joa)
    mean_sigo = np.array(mean_sigo)

    if scaleby != "null":
        mean_scale = np.array(mean_scale)
        mean_ombg = mean_ombg / mean_scale
        rms_ombg = rms_ombg / mean_scale
        mean_oman = mean_oman / mean_scale
        rms_oman = rms_oman / mean_scale
        mean_sigo = mean_sigo / mean_scale

    mean_job = mean_job / sum_nobs
    mean_joa = mean_joa / sum_nobs

    return (
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
