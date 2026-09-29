#Module containing functions specific to tar reading functionality
import os
import tarfile
import tempfile
from pathlib import Path
from dataclasses import replace
from typing import List, Optional, Tuple, Any

from ..loading.iodareader import IODAReader
from ..processing.masking import fill_val_mask, qc_pass_mask, qc_fail_mask
from ..processing.filtering import apply_filter
from ..processing.derived import calc_derived
from .. processing.binning import create_bins
from ..stats.calc_stats import calculate_stats
from ..stats.aggregate import str_to_datetime


#Function to find a specific .nc4 file inside a tarball of a specified instrument(e.g. atms_n20)
def find_instrument_member(tar: tarfile.TarFile, instrument: str) -> Optional[tarfile.TarInfo]:
    members = tar.getmembers()
    #Loop over all members of a tarball
    for member in members:
        base = os.path.basename(member.name)
        
        if base.startswith(instrument) and base.endswith("nc4"):
            return member
         
    return None
    ...

def extract_member_path(tar: tarfile.TarFile, member: tarfile.TarInfo, dest_dir: str) -> str:
     # Reject non-regular files (symlinks, hardlinks, devices, dirs-as-files).
    if not member.isfile():
        raise ValueError(f"Refusing to extract non-regular member: {member.name!r}")

    member_path = os.path.join(dest_dir, member.name)

    # Extract just this member. On Python 3.12+, filter='data' adds another
    # layer of protection; guard so it still works on older versions.
    try:
        tar.extract(member, path=dest_dir, filter="data")  # py3.12+
    except TypeError:
        tar.extract(member, path=dest_dir)                 # older Python

    extracted = os.path.abspath(member_path)
    return extracted

    ...


def _process_single_tar(
    tar_file_path: str,
    instrument: str,
    varname: str,
    kx: int,
    starttime: Any,
    endtime: Any
) -> Optional[Tuple[Any, Any, Any, Any]]:
    # Create the reader inside the worker to avoid pickling issues
    reader = IODAReader()
    exp_name = Path(tar_file_path).name.split(".", 1)[0]

    with tarfile.open(tar_file_path, mode="r:*") as tar:
        member = find_instrument_member(tar, instrument)
        if member is None:
            return None

        with tempfile.TemporaryDirectory() as tmp:
            nc_path = extract_member_path(tar, member, tmp)
            data = reader.read(nc_path, varname, kx)
            data = replace(data, exp=exp_name)

            # Masking and filtering
            val_mask = fill_val_mask(data)
            valid_data = apply_filter(data, val_mask)

            # QC masking
            pass_mask = qc_pass_mask(valid_data)
            fail_mask = qc_fail_mask(valid_data)
            pass_data = apply_filter(valid_data, pass_mask)
            fail_data = apply_filter(valid_data, fail_mask)

            # Derived calculations
            pass_data = calc_derived(pass_data)

            # Binning
            pass_data_binned = create_bins(pass_data)
            pass_data_binned = replace(pass_data_binned, ts_range=[starttime, endtime])
            fail_data_binned = create_bins(fail_data)

            # Stats
            pass_stats_binned = calculate_stats(pass_data_binned)

            return (data.datetime, pass_stats_binned, pass_data_binned, fail_data_binned)
