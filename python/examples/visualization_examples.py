#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from obsview import plot_map, plot_histogram, plot_summary
from obsview.models import ODSData

print("Obsview Python Port - Visualization Examples")
np.random.seed(42)
ods = ODSData()
ods.kt = np.random.choice([4, 5, 11, 33, 44], 100).astype(np.int8)
ods.kx = np.random.choice([120, 220, 281], 100).astype(np.int16)
ods.lon = np.random.uniform(-180, 180, 100).astype(np.float32)
ods.lat = np.random.uniform(-90, 90, 100).astype(np.float32)
ods.lev = np.random.uniform(100, 1000, 100).astype(np.float32)
ods.obs = np.random.normal(0, 5, 100).astype(np.float32)
print("Sample ODS: {} observations".format(len(ods)))
fig1, ax1 = plot_map(ods, domain='global')
print("Global map created")
fig2 = plot_summary(ods)
print("Summary plot created")
plt.close('all')
