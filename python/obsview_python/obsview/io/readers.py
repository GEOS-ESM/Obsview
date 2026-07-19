#Module containing functions that read IODA files and check if a file is an ODS file
#TODO: add a function to read an ODS file in here

"""File-opening logic for the two supported input formats.

- .ods  -> plain netCDF4 files (ODS data-assimilation diagnostics)
- .nc4  -> IODA (JEDI) diagnostic files, optionally packed inside a tarball
"""
import tarfile
import tempfile

from netCDF4 import Dataset

#Function for reading IODA file from a tarball or standalone IODA .nc4 file
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
    netCDF4.Dataset object
    """
    #Check if the user wants the function to open the file as a tarball
    if tarball == "none": #If the user does not specify that the file is a tarball, then return a Dataset() object
        return Dataset(filename, "r")
    
    #If the user does not specify the file is not tarball, then the program assumes it is one and continues 
    print(f"Reading from tarball: {tarball}")
    with tarfile.open(tarball, "r:*") as tar: #Open for reading with transparent compression
        try:
            tarinfo = tar.getmember(filename)
        except KeyError:
            raise FileNotFoundError(f"{filename} not found in tarball.")
        #Create a temporary uncompressed file to read from and return a Dataset() object
        with tempfile.NamedTemporaryFile(suffix=".nc4") as tmp_file:
            fileobj = tar.extractfile(tarinfo)
            tmp_file.write(fileobj.read())
            tmp_file.flush()
            return Dataset(tmp_file.name, "r")

#Returns the lowercase version of the right most part of a file 
#Example: "data.Ods" becomes "ods"
def file_extension(filename):
    """Return the lowercase extension (no dot) used to pick a file format."""
    return filename.rsplit(".", 1)[-1].lower()

#Checks if the end of the given file name ends in .ods
def is_ods(filename):
    return file_extension(filename) == "ods"

#Checks if the file extension ends in .nc4 or .tar
def is_ioda(filename):
    return file_extension(filename) == "nc4" or "tar"
