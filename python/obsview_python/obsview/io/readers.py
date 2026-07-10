"""File-opening logic for the two supported input formats.

- .ods  -> plain netCDF4 files (ODS data-assimilation diagnostics)
- .nc4  -> IODA (JEDI) diagnostic files, optionally packed inside a tarball
"""
import tarfile
import tempfile

from netCDF4 import Dataset


def ioda_from_tarball(tarball, filename):
    """Open an IODA .nc4 file, either directly or extracted from a tarball.

    Parameters
    ----------
    tarball : str
        Path to a tarball containing ``filename``, or the string
        ``"none"`` to open ``filename`` directly.
    filename : str
        Path to the .nc4 file (or its member name inside the tarball).

    Returns
    -------
    netCDF4.Dataset
    """
    if tarball == "none":
        return Dataset(filename, "r")

    print(f"reading from tarball: {tarball}")
    with tarfile.open(tarball, "r:*") as tar:
        try:
            tarinfo = tar.getmember(filename)
        except KeyError:
            raise FileNotFoundError(f"{filename} not found in tarball.")

        with tempfile.NamedTemporaryFile(suffix=".nc4") as tmp_file:
            fileobj = tar.extractfile(tarinfo)
            tmp_file.write(fileobj.read())
            tmp_file.flush()
            return Dataset(tmp_file.name, "r")


def file_extension(filename):
    """Return the lowercase extension (no dot) used to pick a file format."""
    return filename.rsplit(".", 1)[-1].lower()


def is_ods(filename):
    return file_extension(filename) == "ods"
