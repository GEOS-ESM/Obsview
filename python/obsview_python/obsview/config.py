"""Lookup tables that used to live as if/elif chains inside main().

Centralizing these here means adding a new observation type is a data
change (add a dict entry) rather than a code change (add another elif
branch in three different places).
"""

# Default pressure binning range: [bottom, top] in hPa.
DEFAULT_LEVLIM = [1000.0, 0.1]

# obtype -> varname (used when --var is left as 'auto')
OBTYPE_TO_VARNAME = {
    "mls55_aura": "ozoneProfile",
    "gps": "bendingAngle",
    "aircraft_tsen": "airTemperature",
    "temperature": "airTemperature",
    "sondes_tv": "virtualTemperature",
    "sondes_tsen": "airTemperature",
    "sondes_q": "specificHumidity",
    "sondes_u": "windEastward",
    "aircraft_u": "windEastward",
    "sondes_v": "windNorthward",
    "aircraft_v": "windNorthward",
    "saberT": "airTemperature",
    "radiance": "brightnessTemperature",
    "aero": "aerosolOpticalDepth",
}

# obtypes that are binned/plotted by channel number rather than pressure
RADIANCE_OBTYPES = {"radiance"}

# obtype-specific overrides of the default pressure binning range
LEVLIM_OVERRIDES = {
    "radiance": [1, 616],
}

# varname -> kt (JEDI/GSI "kind of table" observation-type code)
VARNAME_TO_KT = {
    "bendingAngle": 89,
    "windEastward": 4,
    "windNorthward": 5,
    "specificHumidity": 11,
    "virtualTemperature": 44,
    "airTemperature": 44,
    "ozoneProfile": 87,
    "brightnessTemperature": 40,
}
DEFAULT_KT = 9999

# --scale CLI value -> internal "scaleby" group name
SCALE_MAP = {
    "obs": "ObsValue",
    "hofx0": "hofx0",
}
DEFAULT_SCALEBY = "null"


class RunConfig:
    """Resolved configuration for a single run, derived from CLI args."""

    def __init__(self, obtype, var, scale):
        self.obtype = obtype
        if var == "auto":
            if obtype not in OBTYPE_TO_VARNAME:
                raise ValueError(
                    f"--obtype {obtype!r} has no default variable; "
                    f"pass --var explicitly (known obtypes: "
                    f"{sorted(OBTYPE_TO_VARNAME)})"
                )
            self.varname = OBTYPE_TO_VARNAME[obtype]
        else:
            self.varname = var
        self.radiance = obtype in RADIANCE_OBTYPES
        self.levlim = list(LEVLIM_OVERRIDES.get(obtype, DEFAULT_LEVLIM))
        self.kt = VARNAME_TO_KT.get(self.varname, DEFAULT_KT)
        self.scaleby = SCALE_MAP.get(scale, DEFAULT_SCALEBY)

    def __repr__(self):
        return (
            f"RunConfig(obtype={self.obtype!r}, varname={self.varname!r}, "
            f"kt={self.kt}, radiance={self.radiance}, levlim={self.levlim}, "
            f"scaleby={self.scaleby!r})"
        )


def resolve(obtype, var, scale):
    """Build a RunConfig from raw CLI argument values."""
    return RunConfig(obtype, var, scale)
