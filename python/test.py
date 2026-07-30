import os
import re
import glob
import tarfile
import tempfile
from typing import Optional

tar_path = "/Users/ltrayano/Desktop/Obsview/Obsview/python/data/IODA files"
instrument_name = "atms_n20"


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


tar_file_paths = sorted(glob.glob(os.path.join(tar_path, "*.tar")))

for tar_file_path in tar_file_paths:
    with tarfile.open(tar_file_path, mode = "r:*") as tar:
        
        member = find_instrument_member(tar, instrument_name)
        print(f"Member is: {member}")
        if member == None:      #If instrument file is missing...
            continue
        with tempfile.TemporaryDirectory() as tmp:
            nc_path = extract_member_path(tar, member, tmp)
            print(f"Temp nc4 path is: {nc_path}")
            ...
        ...