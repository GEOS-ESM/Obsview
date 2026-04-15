"""Visualization module for ODS data."""
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Tuple
import logging
from .models import ODSData

logger = logging.getLogger(__name__)

def plot_map(ods, domain='global', color_by=None, show_grid=True, figsize=None):
    if len(ods) == 0:
        raise ValueError("Cannot plot empty ODS")
    if figsize is None:
        figsize = (12, 8)
    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(ods.lon, ods.lat, s=50, alpha=0.6, c='blue', edgecolors='k', linewidth=0.5)
    if domain == 'global':
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
    elif domain == 'north_america':
        ax.set_xlim(-172, -52)
        ax.set_ylim(16, 72)
    elif domain == 'europe':
        ax.set_xlim(-13, 25)
        ax.set_ylim(35, 72)
    if show_grid:
        ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Observation Locations - {}'.format(domain))
    ax.set_aspect('equal')
    plt.tight_layout()
    return fig, ax

def plot_histogram(ods, attribute, ax=None, bins=20):
    if ax is None:
        fig, ax = plt.subplots()
    data = getattr(ods, attribute)
    valid = np.isfinite(data)
    if not np.any(valid):
        logger.warning("All NaN for {}".format(attribute))
        return ax
    ax.hist(data[valid], bins=bins, alpha=0.7, edgecolor='black')
    ax.set_xlabel(attribute)
    ax.set_ylabel('Count')
    ax.set_title('Distribution of {}'.format(attribute))
    return ax

def plot_summary(ods, figsize=None):
    if figsize is None:
        figsize = (16, 12)
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.scatter(ods.lon, ods.lat, s=20, alpha=0.5)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    ax1.set_title('Observation Locations')
    ax1.grid(True, alpha=0.3)
    ax2 = fig.add_subplot(2, 2, 2)
    plot_histogram(ods, 'kt', ax=ax2)
    ax3 = fig.add_subplot(2, 2, 3)
    plot_histogram(ods, 'kx', ax=ax3)
    ax4 = fig.add_subplot(2, 2, 4)
    plot_histogram(ods, 'lev', ax=ax4, bins=30)
    plt.suptitle('ODS Summary Plot - {} Observations'.format(len(ods)), fontsize=14)
    plt.tight_layout()
    return fig
