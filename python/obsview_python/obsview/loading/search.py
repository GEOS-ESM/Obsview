#Module containing functions to search for ODS and IODA(tar) files 
import os
import glob
from typing import List
from datetime import datetime, timedelta


#Create list of all datetimes based on start and end time
def _synoptic_times(start: datetime, end: datetime) -> List[datetime]:
    times = []
    cursor = start
    step = timedelta(hours=6)
    while cursor <= end:
        times.append(cursor)
        cursor += step
    return times

#Shift time back by 3 hours if the file type is IODA
def _apply_time_convention(dt: datetime, file_type: str) -> datetime:
    if file_type == "ioda":
        return dt - timedelta(hours=3) 
    return dt

#Resolve ODS directory from templates for single synoptic time, 
#find specific file matching instrument and time
def _ods_candidate_file(base_path: str, expid: str, instrument: str,
                         dir_template: str, file_pattern: str,
                         file_time_format: str, dt: datetime) -> List[str]:
   
    # example ODS directory: base/Y2026/M01/D01/H00
    directory = dir_template.format(
        base=base_path,
        expid=expid,
        year=f"{dt.year:04d}",
        month=f"{dt.month:02d}",
        day=f"{dt.day:02d}",
        hour=f"{dt.hour:02d}",
    )
    if not os.path.isdir(directory):
        return []   # this cycle simply has no data; skip silently

    # Filename token for this datetime, e.g. '20260101_00z'.
    time_token = dt.strftime(file_time_format)

    pattern = file_pattern.format(expid=expid, instrument=instrument, file_time_format = time_token)
    
    return sorted(glob.glob(os.path.join(directory, pattern)))


def _ioda_candidate_tar(base_path: str, expid: str, 
                         dir_template: str,
                         file_time_format: str, dt: datetime) -> List[str]:
   
    # example ODS directory: base/Y2026/M01/D01/H00
    directory = dir_template.format(
        base=base_path,
        expid=expid,
        year=f"{dt.year:04d}",
        month=f"{dt.month:02d}")
    
    if not os.path.isdir(directory):
        return []   # this cycle simply has no data; skip silently

    # Filename token for this datetime, e.g. '20260101_00z'.
    time_token = dt.strftime(file_time_format)
    tar_pattern = "{expid}.jedi_hofx.{file_time_format}.tar"
    pattern = tar_pattern.format(expid=expid, file_time_format = time_token)
    
    return sorted(glob.glob(os.path.join(directory, pattern)))













def discover_ods_files(base_path: str, expid: str, instrument: str,
                       dir_template: str, file_pattern: str,
                       file_time_format: str,
                       start: datetime, end: datetime) -> List[str]:

    files = []
    for dt in _synoptic_times(start, end):
        files.extend(_ods_candidate_file(
            base_path, expid, instrument,
            dir_template, file_pattern, file_time_format, dt,
        ))
    return sorted(files)

def discover_ioda_files(base_path: str, expid: str, instrument: str,
                       dir_template: str, file_pattern: str,
                       file_time_format: str,
                       start: datetime, end: datetime) -> List[str]:

    files = []
    for dt in _synoptic_times(start, end):
        dt = _apply_time_convention(dt,"ioda")
        files.extend(_ioda_candidate_tar(base_path,expid,dir_template,file_time_format,dt))
    return sorted(files)




# start = datetime(year = 2026, month=1,day=1,hour=00)
# end = datetime(year = 2026, month=1,day=1,hour=18)
# base_path = "/Users/ltrayano/Desktop/Obsview/Obsview/python/data/IODA files"
# expid = "j54rp1"
# dir_template = "{base}/{expid}/jedi/obs/Y{year}/M{month}"
# file_time_format = "%Y%m%d_%Hz"



# times = _synoptic_times(start,end)

# files = []
# for time in times:
#     time = _apply_time_convention(time,"ioda")
#     files.extend(_ioda_candidate_tar(base_path,expid,dir_template,file_time_format,time))

# print(files)