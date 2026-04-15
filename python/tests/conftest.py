"""Pytest configuration for obsview tests"""
import pytest
import numpy as np
from obsview.models import ODSData

@pytest.fixture
def minimal_ods():
    ods = ODSData()
    ods.kt = np.array([4, 5, 33], dtype=np.int8)
    ods.kx = np.array([220, 220, 181], dtype=np.int16)
    ods.lon = np.array([0, 90, -90], dtype=np.float32)
    ods.lat = np.array([0, 45, -30], dtype=np.float32)
    ods.lev = np.array([500, 700, 1000], dtype=np.float32)
    ods.obs = np.array([10.0, 5.0, 1013.0], dtype=np.float32)
    return ods

@pytest.fixture
def large_ods():
    np.random.seed(42)
    nobs = 1000
    ods = ODSData()
    ods.kt = np.random.choice([4, 5, 11, 33, 44], nobs).astype(np.int8)
    ods.kx = np.random.choice([120, 220, 281, 181], nobs).astype(np.int16)
    ods.lon = np.random.uniform(-180, 180, nobs).astype(np.float32)
    ods.lat = np.random.uniform(-90, 90, nobs).astype(np.float32)
    ods.lev = np.random.uniform(100, 1000, nobs).astype(np.float32)
    ods.obs = np.random.normal(0, 5, nobs).astype(np.float32)
    return ods
