#Module containing the function to combine U and V wind components
import numpy as np
from dataclasses import replace
from.filtering import apply_filter
from .masking import uv_match_masks
from ..loading.observationdata import ObservationData


def combine_uv_components(u: ObservationData, v: ObservationData) -> ObservationData:
    #New arrays with matching observations:
    iu, iv = uv_match_masks(u,v)
    u,v = apply_filter(u,iu), apply_filter(v,iv)

    #Background/analysis components
    ub, vb = u.obs - u.omb, v.obs - v.omb
    ua, va = u.obs - u.oma, v.obs - v.oma

    #Obs magnitude
    mag_obs = np.hypot(u.obs,v.obs)

    #Residual vectors
    omb = mag_obs - np.hypot(ub, vb)
    oma = mag_obs - np.hypot(ua, va)

    amb = np.hypot(u.amb, v.amb)                         # |analysis - background|
    sigo = np.sqrt(0.5 * (u.sigo**2 + v.sigo**2))        # == sigo when U,V share it

    # Vector passes QC only if BOTH components pass; otherwise carry the failure.
    both_pass = (u.qc == 0) & (v.qc == 0)
    qc = np.where(both_pass, 0, np.maximum(u.qc, v.qc))

    # Estimated-error cross terms must be dot products, not |omb|*|oma|.
    esigo = u.omb * u.oma + v.omb * v.oma
    esigb = u.omb * u.amb + v.omb * v.amb

    # All scalar metadata (datetime, lev_type, file_type, all_lev, fill_values,
    # exp) and the per-obs lat/lon/lev/kx come from the matched U records.
    return replace(u, obs=mag_obs, omb=omb, oma=oma, amb=amb, sigo=sigo,
                   qc=qc, esigo=esigo, esigb=esigb)
    ...