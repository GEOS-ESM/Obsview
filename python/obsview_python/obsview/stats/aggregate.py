#Module for aggregating data from TimeSeriesData objects
import re
import numpy as np
from datetime import datetime, timedelta, timezone
from typing import List

from ..loading.observationdata import ObservationData
from ..loading.timeseriesdata import TimeSeriesData
from .statisticsdata import StatisticsData
from ..processing.binning import BinnedData

#TODO: Put this in utils.py
#Create datetime object from string formatted like Matlab version:
# expected string format: 'YYYYMMDDHH', example: '2026010100'
def str_to_datetime(datetime_str:str) -> datetime:
    m = re.search(r"(\d{8})(\d{2})", datetime_str, flags=re.IGNORECASE)
    date_str, hour_str = m.group(1), m.group(2)
    dt = datetime.strptime(date_str + hour_str, "%Y%m%d%H")
    return dt.replace(tzinfo=timezone.utc)  
    ...

#Create list of datetime objects from user specified start and end datetime strings
def create_time_range(start_time:str, end_time:str) -> List[datetime]:
    starttime = str_to_datetime(start_time)
    endtime = str_to_datetime(end_time)

    # Calculate total 6-hour intervals
    steps = int((endtime - starttime).total_seconds() / 3600) // 6

    # Generate list of datetime objects
    synoptic_datetimes = [starttime + timedelta(hours=i * 6) for i in range(steps + 1)]
    return synoptic_datetimes
    ...

#Function to compute time averaged statistics from user specified time range (starttime, endtime)
def aggregate_stats(ts: TimeSeriesData, start_time: str, end_time: str) -> StatisticsData:
    lev_count = np.size(ts.pass_stats[0].nobs)
    
    time_range = create_time_range(start_time, end_time)
    
    dt_count = 0
    nobs = np.zeros(lev_count)
    mean_omb = np.zeros(lev_count)
    rms_omb = np.zeros(lev_count)
    mean_oma = np.zeros(lev_count)
    rms_oma = np.zeros(lev_count)
    mean_job = np.zeros(lev_count)
    mean_joa = np.zeros(lev_count)
    mean_sigo = np.zeros(lev_count)
    mean_esigo = np.zeros(lev_count)
    mean_esigb = np.zeros(lev_count)

    for datetime in time_range:     #Loop over all dts in time range
        dt_count += 1
        if datetime in ts.datetimes:    #If specific datetime from range is equal to one of the timeseries datetimes...
            nobs += ts.pass_stats[ts.datetimes.index(datetime)].nobs
            mean_omb =+ ts.pass_stats[ts.datetimes.index(datetime)].mean_omb
            rms_omb += ts.pass_stats[ts.datetimes.index(datetime)].rms_omb
            mean_oma += ts.pass_stats[ts.datetimes.index(datetime)].mean_oma
            rms_oma += ts.pass_stats[ts.datetimes.index(datetime)].rms_oma
            mean_job += ts.pass_stats[ts.datetimes.index(datetime)].mean_job
            mean_joa += ts.pass_stats[ts.datetimes.index(datetime)].mean_joa
            mean_sigo += ts.pass_stats[ts.datetimes.index(datetime)].mean_sigo
            mean_esigo += ts.pass_stats[ts.datetimes.index(datetime)].mean_esigo
            mean_esigb += ts.pass_stats[ts.datetimes.index(datetime)].mean_esigb

    #Calculate time averaged variables
    avg_nobs = nobs/dt_count
    avg_mean_omb = mean_omb/dt_count
    avg_rms_omb = rms_omb/dt_count
    avg_mean_oma = mean_oma/dt_count
    avg_rms_oma = rms_oma/dt_count
    avg_mean_job = mean_job/dt_count
    avg_mean_joa = mean_joa/dt_count
    avg_mean_sigo = mean_sigo/dt_count
    avg_mean_esigo = mean_esigo/dt_count
    avg_mean_esigb = mean_esigb/dt_count
    
    #Assign time averaged variables to StatisticsData attibutes

    obj = StatisticsData(
        nobs = avg_nobs,
        mean_omb = avg_mean_omb,
        rms_omb = avg_rms_omb,
        mean_oma = avg_mean_oma,
        rms_oma = avg_rms_oma,
        mean_job = avg_mean_job,
        mean_joa = avg_mean_joa,
        mean_sigo = avg_mean_sigo,
        mean_esigo = avg_mean_esigo,
        mean_esigb = avg_mean_esigb
    )

    return obj


def aggregate_pass_binned(ts: TimeSeriesData, start_time:str, end_time: str) -> BinnedData:
    lev_count = np.size(ts.pass_stats[0].nobs)
    time_range = create_time_range(start_time, end_time)
    
    
    #For data object, keep:
    lev_type = ts.pass_data[0].data.lev_type
    kx = ts.pass_data[0].data.kx
    kt = ts.pass_data[0].data.kt
    file_type = ts.pass_data[0].data.file_type

    #Pass only the data that gets used for making statistics plot
    data = ObservationData(
        lev_type = lev_type,
        kx = kx,
        kt = kt,
        file_type = file_type
    )
    
    #For binned data object, aggregate:
    #bin indices
    dt_count = 0
    nobs = np.zeros(lev_count)
    for datetime in time_range:     #Loop over all dts in time range
        dt_count += 1
        if datetime in ts.datetimes:
            nobs += np.bincount(ts.pass_data[ts.datetimes.index(datetime)].bin_indices)
             
    avg_nobs = nobs/dt_count
    bin_centers = ts.pass_data[0].bin_centers
    bin_labels = ts.pass_data[0].bin_labels
    bin_heights = ts.pass_data[0].bin_heights
    is_ts = True
    starttime = str_to_datetime(start_time)
    endtime = str_to_datetime(end_time)
    ts_range = [starttime, endtime]

    obj = BinnedData(
        data = data,
        bin_centers = bin_centers,
        bin_labels = bin_labels,
        bin_heights = bin_heights,
        nobs = avg_nobs,
        is_ts = is_ts,
        ts_range = ts_range
    )
    return obj
    

def aggregate_fail_binned(ts: TimeSeriesData, start_time:str, end_time: str) -> BinnedData:
    lev_count = np.size(ts.pass_stats[0].nobs)
    time_range = create_time_range(start_time, end_time)
    
    
    #For data object, keep:
    lev_type = ts.pass_data[0].data.lev_type
    kx = ts.pass_data[0].data.kx
    kt = ts.pass_data[0].data.kt

    #Pass only the data that gets used for making statistics plot
    data = ObservationData(
        lev_type = lev_type,
        kx = kx,
        kt = kt
    )
    
    #For binned data object, aggregate:
    #bin indices
    dt_count = 0
    nobs = np.zeros(lev_count)
    for datetime in time_range:     #Loop over all dts in time range
        dt_count += 1
        if datetime in ts.datetimes:
            nobs += np.bincount(ts.fail_data[ts.datetimes.index(datetime)].bin_indices)
             
    avg_nobs = nobs/dt_count
    bin_centers = ts.fail_data[0].bin_centers
    bin_labels = ts.fail_data[0].bin_labels
    bin_heights = ts.fail_data[0].bin_heights

    obj = BinnedData(
        data = data,
        bin_centers = bin_centers,
        bin_labels = bin_labels,
        bin_heights = bin_heights,
        nobs = avg_nobs
    )
    return obj
    ...
       