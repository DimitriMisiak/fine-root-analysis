#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 11:52:40 2020

@author: misiak
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from representation_functions import nodecor_crosstalk_correction, ax_hist
from data_analysis import (
    glitch_time_cut,
    heat_chi2_cut,
    ion_chi2_cut,
    offset_ion_cut,
    quality_cut
)
    

fine_data_path = '/home/misiak/Analysis/neutron_background/data_analysis.h5'

source = 'Background'
DF = pd.read_hdf(
    fine_data_path,
    where='source = "{}"'.format(source))

#%%
plt.close('all')

fine_cut = DF['glitch_time_cut'] & DF['quality_cut'] #& DF['bulk_cut']

bins = np.arange(-10, 10, 0.1)


fig, ax = plt.subplots()
ax_hist(
        ax,
        bins,
        DF[fine_cut]['energy_nodecor_ion_conservation'],
        'nodecor'
)


# ax_hist(
#         ax,
#         bins,
#         DF[fine_cut]['energy_ion_conservation'],
#         'decor',
#         color='coral',
# )

energy_heat = DF[fine_cut]['energy_heat']
# energy_heat = DF[fine_cut]['recoil_energy_total']
ion_conservation = DF[fine_cut]['energy_nodecor_ion_conservation']

from scipy.stats import binned_statistic

bins = np.arange(0, 1000, 100)
bins_width = bins[1] - bins[0]
bins_array = bins[:-1] + (bins_width) / 2

bins = np.logspace(np.log10(0.2), np.log10(500), 25)
bins_width = (bins[1:] - bins[:-1])
bins_array = bins[:-1]

def std_like(x_array):
    inf, sup = np.quantile(x_array, [0.16, 0.84])
    return (sup - inf)/2

stats, bin_edges, bin_number = binned_statistic(
    energy_heat,
    ion_conservation,
    statistic=std_like,
    bins=bins
)

from scipy.optimize import curve_fit

def std_model(x_array, a, b):
    return a + b * x_array

popt, pcov = curve_fit(std_model, bins_array, stats, sigma=0.1*stats)

fig, ax = plt.subplots()
ax.plot(bins_array, stats, ls='steps-mid')

x_array = np.linspace(0, 500, int(1e4))
mod = std_model(x_array, *popt)
ax.plot(x_array, mod)

means, bin_edges, bin_number = binned_statistic(
    energy_heat,
    ion_conservation,
    statistic=(lambda x: np.quantile(x, [0.5,])),
    bins=bins,
)

fig, ax = plt.subplots()
ax.plot(energy_heat, ion_conservation, ls='none', marker='.', alpha=0.1)

ax.errorbar(
    bins_array,
    means,
    yerr=stats,
    marker='o',
    color='k',
    zorder=10
)
