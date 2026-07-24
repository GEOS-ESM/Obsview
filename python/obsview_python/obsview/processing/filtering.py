#Module to contain functions that keep/hide data based on 1-D mask arrays
import numpy as np
from dataclasses import fields, replace
from ..loading.observationdata import ObservationData


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