# Configuration file containing settings category dataclass objects
# and yaml file loading/validating functions
import yaml
from dataclasses import dataclass, field
from typing import List, Optional

# ============================================================
# Runtime settings
# ============================================================
@dataclass
class RuntimeSettings:
    parallel: bool = True
    verbose: bool = True

# ------------------------------------------------------------
# Data source templates
# ------------------------------------------------------------
@dataclass
class SourceTemplate:
    dir_template: str
    file_time_format: str
    file_pattern: str


@dataclass
class DataTemplate:
    ioda: SourceTemplate
    ods: SourceTemplate

# ------------------------------------------------------------
# Experiment locations
# ------------------------------------------------------------
@dataclass
class Experiment:
    id: str
    file_type: str
    base_path: str
    label: str

# ------------------------------------------------------------
# Analysis selection
# ------------------------------------------------------------
@dataclass
class TimeRange:
    start: str
    end: str

@dataclass
class Target:
    instrument: str
    kx: int
    obtype: str
    varname: str
    kt: int
    level: int
    
    
    

@dataclass
class Selection:
    time_range: TimeRange
    targets: List[Target]
    domain: str

# ------------------------------------------------------------
# Plot settings
# ------------------------------------------------------------
@dataclass
class StatisticsSettings:
    enabled: bool
    time_averaged: bool


@dataclass
class ComparisonSettings:
    stat_type: str
    ratio_bounds: float
    enabled: bool = True

@dataclass
class CoverageMapSettings:
    enabled: bool
    map_level: int
    projection: str

@dataclass
class TimeSeriesSettings:
    enabled: bool
    date_tick_interval_days: int


@dataclass
class PlotSettings:
    statistics: StatisticsSettings
    comparison: ComparisonSettings
    coverage_map: CoverageMapSettings
    time_series: TimeSeriesSettings

# ------------------------------------------------------------
# Output settings
# ------------------------------------------------------------
@dataclass
class Output:
    mode: str
    directory: str
    format: str
    dpi: int
    filename_template: str




#Overarching config dataclass
@dataclass
class Config:
    run: RuntimeSettings
    data: DataTemplate
    experiments: List[Experiment]
    selection: Selection
    plots: PlotSettings                      # or a typed PlotSettings
    output: dict                     # or a typed OutputSettings











#Load config file contents into Config object
def load_config(path: str) -> Config:
    with open(path, "r") as f:
        raw = yaml.safe_load(f)          # SAFE loader only

    _validate(raw)                        # fail loud before building objects

    #run
    run = RuntimeSettings(**raw.get("run", {}))
    #data
    data = DataTemplate(
        ioda=SourceTemplate(**raw["data"]["ioda"]),
        ods=SourceTemplate(**raw["data"]["ods"]),
    )
    #experiements
    experiments = [Experiment(**e) for e in raw["experiments"]]
    #selection
    sel_raw = raw["selection"]
    selection = Selection(
        time_range=TimeRange(**sel_raw["time_range"]),
        targets=[Target(**t) for t in sel_raw["targets"]],
        domain=sel_raw.get("domain"),     # optional
    )
    #plots
    plots = PlotSettings(
        statistics= StatisticsSettings(**raw["plots"]["statistics"]),
        comparison = ComparisonSettings(**raw["plots"]["comparison"]),
        coverage_map = CoverageMapSettings(**raw["plots"]["coverage_map"]),
        time_series = TimeSeriesSettings(**raw["plots"]["time_series"])
    )
    #output
    output = Output(**raw["output"])

    return Config(
        run=run,
        data=data,
        experiments=experiments,
        selection=selection,
        plots=plots,
        output=output,
    )

def _validate(raw: dict) -> None:
    """Check required keys / allowed values; raise with a clear message."""

    if not raw.get("experiments"):
        raise ValueError("At least one experiment must be defined")
    if not raw["selection"].get("targets"):
        raise ValueError("selection.targets must contain at least one target")
    # ... more checks ...