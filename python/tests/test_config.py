"""Tests for configuration module"""
import pytest
from obsview.config import dconfig, ODSConfig

def test_dconfig_single_parameter():
    kts = dconfig('KTS')
    assert kts is not None

def test_dconfig_multiple_parameters():
    kts, kxs = dconfig('KTS', 'KXS')
    assert kts is not None and kxs is not None

def test_dconfig_invalid_parameter():
    with pytest.raises(ValueError):
        dconfig('INVALID_PARAM')

def test_ods_config_initialization():
    config = ODSConfig()
    assert config.data_types is not None
