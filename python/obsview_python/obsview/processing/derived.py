#Module for calculating variable derived from ODS or IODA file variables
from dataclasses import replace
from ..loading.observationdata import ObservationData

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
