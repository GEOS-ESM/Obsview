"""Tests for data loading module"""
import pytest
import numpy as np
from obsview.models import ODSData

def test_odsdata_creation():
    ods = ODSData()
    assert len(ods) == 0

def test_odsdata_with_data():
    ods = ODSData()
    ods.kt = np.array([4, 5], dtype=np.int8)
    ods.lon = np.array([0, 90], dtype=np.float32)
    assert len(ods) == 2

def test_odsdata_subset():
    ods = ODSData()
    ods.kt = np.array([4, 5, 4, 33], dtype=np.int8)
    mask = ods.kt == 4
    ods_subset = ods.subset(mask)
    assert len(ods_subset) == 2
