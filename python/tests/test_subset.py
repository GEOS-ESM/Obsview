"""Tests for subsetting module"""
import pytest
import numpy as np
from obsview.subset import odssubset, odsclean
from obsview.models import ODSData

def test_odssubset_boolean_mask():
    ods = ODSData()
    ods.kt = np.array([4, 5, 4, 33], dtype=np.int8)
    mask = ods.kt == 4
    ods_sub = odssubset(ods, mask)
    assert len(ods_sub) == 2

def test_odsclean():
    ods = ODSData()
    ods.kt = np.array([4, 5, 4], dtype=np.int8)
    ods.qcx = np.array([0, 1, 0], dtype=np.int8)
    ods_clean = odsclean(ods)
    assert len(ods_clean) == 2
