#!/usr/bin/env python3
from obsview import odsload
from obsview.config import dconfig
import numpy as np
from obsview.models import ODSData

print("Obsview Python Port - Basic Loading Example")
kts = dconfig('KTS')
print("Data types configured: {}".format(len(kts)))
ods = ODSData()
ods.filename = 'sample.ods'
ods.kt = np.array([4, 5, 33], dtype=np.int8)
ods.kx = np.array([220, 220, 181], dtype=np.int16)
ods.lon = np.array([0, 90, -90], dtype=np.float32)
ods.lat = np.array([0, 45, -30], dtype=np.float32)
print("Loaded: {}".format(ods))
print("Example complete!")
