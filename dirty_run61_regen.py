#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:22:11 2019

@author: misiak
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from spec_classes import Analysis_red
from representation import (
        temporal_plot, plot_chi2_vs_energy,
        histogram_adu, histogram_ev, ion_vs_ion,
        virtual_vs_virtual_ev, optimization_info,
)
from model_spectrum import fid_mixture, double_norm
from plot_addon import LegendTitle, custom_autoscale, ax_hist, plot_ion_vs_ion
from stats_addon import cdf_calc, custom_bin_edges

run_dir = '/home/misiak/Data/data_run61'
              

detector = 'REDN1'


# first command
plt.close('all')
plt.rcParams['text.usetex']=True

n_sigma = 5

### for tg12l003
stream_16 = 'tk16l000'  
ana_16 = Analysis_red(
        stream_16,
        detector=detector,
        run_dir=run_dir,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=n_sigma,
        override_mu=(861, -54, 54),
        override_sigma=(40, 10, 9)
)

stream_15 = 'tk15l005'  
ana_15 = Analysis_red(
        stream_15,
        detector=detector,
        run_dir=run_dir,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=n_sigma,
        override_mu=(861, -54, 54),
        override_sigma=(40, 10, 9)
)

#%%

save_dir = '/home/misiak/Analysis/fine_root_analysis/crude_run61'

fpath_15 = '/'.join((save_dir, stream_15+'.npy'))
fpath_16 = '/'.join((save_dir, stream_16+'.npy'))

crude_15 = np.load(fpath_15)
crude_16 = np.load(fpath_16)

time_15 = ana_15.all.trig.time[ana_15.all.trig.cut.quality][crude_15]
time_16 = ana_16.all.trig.time[ana_16.all.trig.cut.quality][crude_16]
time_tot = np.concatenate([time_15, time_16])


fpath_15 = '/'.join((save_dir, stream_15+'_selection_bulk.npy'))

bulk_15 = np.load(fpath_15)
bulk_16 = np.load(fpath_16)
time_bulk = np.concatenate([time_15[bulk_15], time_16[bulk_16]])


fpath_15 = '/'.join((save_dir, stream_15+'_selection_frontier.npy'))

front_15 = np.load(fpath_15)
front_16 = np.load(fpath_16)
time_frontier = np.concatenate([time_15[front_15], time_16[front_16]])

#time_10kev = ana.all.trig.time[ana.all.trig.cut.quality][crude_ind][bulk_ind]
#energy = energy_raw[crude_ind][bulk_ind, 2:].T

#time_tot = np.concatenate([time_tot_tk15, time_tot_tk16])
#
#time_frontier = np.concatenate([time_front_tk15, time_front_tk16])
#
#time_10kev = np.concatenate([time_10kev_tk15, time_10kev_tk16])
#energy = np.concatenate([energy_tk15, energy_tk16], axis=1)

fig, axes = plt.subplots(ncols=4, num='ion vs time', figsize=(15,6))
axes[0].grid()


n, bins, patches = axes[1].hist(time_bulk, bins=10)

n, bins2, patches = axes[2].hist(time_tot, bins=10)

n, bins3, patches = axes[3].hist(time_frontier, bins=10)

#funk = lambda t,a,b: a + b*t
#from scipy.optimize import curve_fit
#
#for j, ion in enumerate(energy):
#    if j==1 or j==3:
##        axes[0].axhline(np.mean(ion), ls='--', color='k', zorder=10)
#        y_data = (ion-np.mean(ion))/np.mean(ion)
#        axes[0].plot(time_10kev, y_data, ls='none', marker='.')
#        popt, pcov = curve_fit(funk, time_10kev, y_data, p0=[np.mean(ion), 0])
#        axes[0].plot(time_10kev, funk(time_10kev, *popt), 'k')
#        
ax0 = ax_hist(axes[1], bins, time_bulk, 'Timestamp')[0]
axes[1].set_yscale('log')
ax0.plot(time_bulk[[0, -1]],[0, 1], color='k')

ax1 = ax_hist(axes[2], bins2, time_tot, 'Timestamp')[0]
axes[2].set_yscale('log')
ax1.plot(time_tot[[0, -1]],[0, 1], color='k')

ax2 = ax_hist(axes[3], bins3, time_frontier, 'Timestamp')[0]
axes[3].set_yscale('log')
ax2.plot(time_frontier[[0, -1]],[0, 1], color='k')
