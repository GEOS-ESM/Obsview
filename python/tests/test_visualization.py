"""Tests for visualization module"""
import pytest
import matplotlib.pyplot as plt
from obsview.visualization import plot_map, plot_histogram, plot_summary
from obsview.models import ODSData

def test_plot_map_basic(minimal_ods):
    fig, ax = plot_map(minimal_ods)
    assert fig is not None
    plt.close(fig)

def test_plot_map_empty_raises():
    ods = ODSData()
    with pytest.raises(ValueError):
        plot_map(ods)

def test_plot_histogram_kt(minimal_ods):
    fig, ax = plt.subplots()
    plot_histogram(minimal_ods, 'kt', ax=ax)
    plt.close(fig)

def test_plot_summary(minimal_ods):
    fig = plot_summary(minimal_ods)
    assert fig is not None
    plt.close(fig)
