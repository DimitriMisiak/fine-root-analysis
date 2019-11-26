#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
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
stream = 'tk20l001'  
ana_1 = Analysis_red(
        stream,
        detector=detector,
        run_dir=run_dir,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=n_sigma,
        override_mu=(861, -54, 54),
        override_sigma=(40, 10, 9)
)

#stream = 'tk15l005'  
#ana_2 = Analysis_red(
#        stream,
#        detector=detector,
#        run_dir=run_dir,
#        chan_valid=(0, 2, 3, 4, 5),
#        chan_signal=(0, 3, 5),
#        n_sigma=n_sigma,
#        override_mu=(861, -54, 54),
#        override_sigma=(40, 10, 9)
#)

#%%

plt.close('all')

##### USUAL PLOTSSS
ana=ana_1
fig_temp = temporal_plot(ana)
fig_chi2_trig, fig_chi2_noise = plot_chi2_vs_energy(ana)
fig_hist_trig, fig_hist_noise = histogram_adu(ana)


#%%
### CUSTOM ANALYSIS
import red_magic as rmc

# sensitvity estimation a la mano
fig_ion = ion_vs_ion(ana, n_sigma)

axes = fig_ion.get_axes()
ax = axes[3]
ax.set_title('')

line = ax.lines[1]
print(line.get_label())

def stats_funk(indexes):
    x_data = line.get_xdata()[indexes]
    y_data = line.get_ydata()[indexes]
    
    print('X Axis: mean={}    and    std={}'.format(x_data.mean(), x_data.std()))
    print('Y Axis: mean={}    and    std={}'.format(y_data.mean(), y_data.std()))

DS = rmc.Data_Selector(ax, line, proceed_func=stats_funk)

#%%
pos10_dict = {
        'tk18l000':[300, 1, 50.20, -56.44, -52.32, 56.86],
        'tk15l005':[383, 1, 50.56, -57.12, -51.18, 57.63],
        'tk16l000':[371, 1, 48.93, -55.28, -52.44, 56.42],
        'tk18l001':[395, 1, -50.52, 59.07, 51.83, -58.01],
        'tk19l000':[265, 1, 47.55, -49.86, -50.31, 53.43],
        'tk19l001':[233, 1, 39, -47.21, -42.86, 46.63],
        'tk20l000':[562, 1, 49.45, -59.36, -51.95, 59.26],
        'tk20l001':[541, 1, -55, -55, 55, 55],
        'tk20l003':[350, 1, 52.86, -51.2, -53.93, 53.10],
}


# calibration a la mano
energy_raw = ana.all.trig.filt_decor.Energy_OF[ana.all.trig.cut.quality]

energy_calib = energy_raw * 10.37 / np.array(pos10_dict[stream])

heat = energy_calib[:,0]
ion_tot = np.sum(energy_calib[:,2:], axis=1)/2

fig, ax = plt.subplots(num='Tot Ion vs Heat')
ax.plot(heat, ion_tot, label='quality events', ls='none', marker='.')
ax.grid()

save_dir = '/home/misiak/Analysis/fine_root_analysis/crude_run61'
fpath = '/'.join((save_dir, stream+'.npy'))

def crude_funk(indexes):
    np.save(fpath, indexes)
    print('Indexes were saved !')

line = ax.lines[0]
DS_2 = rmc.Data_Selector(ax, line, proceed_func=crude_funk)

#%%

crude_ind = np.load(fpath)

ax.plot(heat[crude_ind], ion_tot[crude_ind], label='quality events', ls='none', marker='.', color='r')

energy_10kev = energy_calib[crude_ind]

heat_adv, _, ionA, ionB, ionC, ionD = energy_10kev.T

fig = plot_ion_vs_ion(ana, energy_10kev, marker=',')

axes = fig.get_axes()

ax = axes[3]
ax.set_title('gotcha')

print('Total number: ', crude_ind.size)

fpath_2 = '/'.join((save_dir, stream+'_selection.npy'))
def count_funk(indexes):
    np.save(fpath_2, indexes)
    print('Number in Selection: ', indexes.size)

line = ax.lines[0]
DS_3 = rmc.Data_Selector(ax, line, proceed_func=count_funk)

#%%
##Checking the regen: ion in function of the time
#bulk_ind = np.load(fpath_2)
#
##time_10kev = ana.all.trig.time[ana.all.trig.cut.quality][crude_ind][bulk_ind]
##energy = energy_raw[crude_ind][bulk_ind, 2:].T
#
#time_tot = np.concatenate([time_tot_tk15, time_tot_tk16])
#
#time_frontier = np.concatenate([time_front_tk15, time_front_tk16])
#
#time_10kev = np.concatenate([time_10kev_tk15, time_10kev_tk16])
#energy = np.concatenate([energy_tk15, energy_tk16], axis=1)
#
#fig, axes = plt.subplots(ncols=4, num='ion vs time', figsize=(15,6))
#axes[0].grid()
#
#
#n, bins, patches = axes[1].hist(time_10kev, bins=10)
#
#n, bins2, patches = axes[2].hist(time_tot, bins=10)
#
#n, bins3, patches = axes[3].hist(time_frontier, bins=10)
#
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
#ax0 = ax_hist(axes[1], bins, time_10kev, 'Timestamp')[0]
#axes[1].set_yscale('log')
#ax0.plot(time_10kev[[0, -1]],[0, 1], color='k')
#
#ax1 = ax_hist(axes[2], bins2, time_tot, 'Timestamp')[0]
#axes[2].set_yscale('log')
#ax1.plot(time_tot[[0, -1]],[0, 1], color='k')
#
#ax2 = ax_hist(axes[3], bins3, time_frontier, 'Timestamp')[0]
#axes[3].set_yscale('log')
#ax2.plot(time_frontier[[0, -1]],[0, 1], color='k')
#
