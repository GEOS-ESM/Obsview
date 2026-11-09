#Module defining the Reference object and respective functions/subclasses
import yaml
from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class Instrument:
    name: str
    title: str
    obtype: str


@dataclass
class Variable:
    name: str
    kt: int
    title: str


@dataclass
class Observation:
    instrument: str
    varname: str
    kx: int
    level_type: str
    levels: str


@dataclass
class Reference:
    instruments: Dict[str,Instrument] = field(default_factory=dict)
    variables: Dict[str,Variable] = field(default_factory=dict)
    observations: List[Observation] = field(default_factory=list)
    level_sets: Dict[str,List[float]] = field(default_factory=dict)

    def get_observation(self, instrument: str, varname: str, kx: Optional[int] = None) -> Optional[Observation]:
        for obs in self.observations:
            if obs.instrument == instrument and obs.varname == varname:
                if kx is None or obs.kx == kx:
                    return obs
        return None
    
    def get_kt_from_varname(self, varname: str) -> int:
         kt  = self.variables.get(varname).kt
         return kt
         ...

    def get_pressure_levels(self, obs: Observation) -> Optional[List[float]]:
        if obs.level_type == 'pressure' and obs.levels:
            if isinstance(obs.levels, str):
                # It's a reference to a level set
                return self.level_sets.get(obs.levels)
            elif isinstance(obs.levels, list):
                # It's already a list of levels
                return obs.levels
        return None
    
    def get_instrument_title(self, instrument: str) -> str:
        inst_title = self.instruments.get(instrument).title
        return inst_title
    
    def get_variable_title(self, variable: str) -> str:
        var_title = self.variables.get(variable).title
        return var_title






def load_reference(path: str) -> Reference:
    with open(path, 'r') as f:
        data = yaml.safe_load(f)

    ref = Reference()

    ref.level_sets = data['level_sets']

    for name, attrs in data['instruments'].items():
            ref.instruments[name] = Instrument(
                name=name,
                title=attrs['title'],
                obtype=attrs['obtype']
            )

    for name, attrs in data['variables'].items():
            ref.variables[name] = Variable(
                name=name,
                kt=attrs['kt'],
                title=attrs['title']
            )

    for obs_data in data['observations']:
            ref.observations.append(Observation(
                instrument=obs_data['instrument'],
                varname=obs_data['varname'],
                kx=obs_data['kx'],
                level_type=obs_data['level_type'],
                levels=obs_data.get('levels')
            )) 

    return ref

