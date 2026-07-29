from dataclasses import dataclass
from typing import List
from ..stats.statisticsdata import StatisticsData
from ..processing.binning import BinnedData

@dataclass
class TimeSeriesData:
    datetimes: List[object]
    pass_stats: List[StatisticsData]
    pass_data: List[BinnedData]
    fail_data: List[BinnedData]
    ...

