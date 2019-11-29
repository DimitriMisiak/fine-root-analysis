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
import red_magic as rmc


# matplotlib configuration
plt.close('all')
plt.rcParams['text.usetex']=True


stream = 'tk25l001'
save_flag = True

# global variables
SAVE_DIR = '/home/misiak/Analysis/fine_root_analysis/crude_run61/{}'.format(stream)
DATA_DIR = '/home/misiak/Data/data_run61'
os.makedirs(SAVE_DIR, exist_ok=True)

# getting the data
ana = Analysis_red(
        stream,
        detector='REDN1',
        run_dir=DATA_DIR,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=5,
        override_mu=(861, -54, 54),
        override_sigma=(40, 10, 9)
)


# =============================================================================
# PLOTS
# =============================================================================

# monitoring vs TIME
fig_temp = temporal_plot(ana)

# QUALITY CUTS
fig_chi2_trig, fig_chi2_noise = plot_chi2_vs_energy(ana)

# resize the plots
for ax in fig_chi2_trig.get_axes():
    ax.set_xlim(10**-2, 10**5)
    ax.set_ylim(10**1, 10**9)

# HEAT CALIBRATION
fig_hist_trig, fig_hist_noise = histogram_adu(ana)

# resize the plots
fig_hist_trig.get_axes()[0].set_xlim(-200, 2000)

if save_flag:
    fig_temp.savefig(SAVE_DIR+'/fig_temp.png')
    fig_chi2_trig.savefig(SAVE_DIR+'/fig_chi2_trig.png')
    fig_chi2_noise.savefig(SAVE_DIR+'/fig_chi2_noise.png')
    fig_hist_trig.savefig(SAVE_DIR+'/fig_hist_trig.png')
    fig_hist_noise.savefig(SAVE_DIR+'/fig_hist_noise.png')   

#%%
# ION CALIBRATION B+D
fig_ion = ion_vs_ion(ana, ana.n_sigma)

# HACK to plot multiples figures using the "ion_vs_ion" function
label = 'ION CALIBRATION B+D'
fig_ion.set_label(label)
fig_ion.canvas.set_window_title(label)

axes = fig_ion.get_axes()

for ax in axes:
    ax.set_xlim(-80, 80)
    ax.set_ylim(-80, 80)

if save_flag:
    fig_ion.savefig(SAVE_DIR+'/fig_ion_all.png')

fig_ion.set_size_inches(15, 10, forward=True)

# selecting the ionB vs ionD plot
ax = axes[3]
# if lost, modify the title to vizualize the ax in which to make the selection
ax.set_title('')

# HACK get the line corresponding to the quality events
line = ax.lines[1]
assert (line.get_label() == 'Quality Events')

def stats_funk(indexes):
    """
    Used for the proceed button in the Data_selector class.
    Just print the mean and std of the selection along the X and Y axis.
    """
    x_data = line.get_xdata()[indexes]
    y_data = line.get_ydata()[indexes]
    print((
            'X Axis: mean={}    and    std={}'
    ).format(x_data.mean(), x_data.std()))
    print((
            'Y Axis: mean={}    and    std={}'
    ).format(y_data.mean(), y_data.std()))

# calling the Data_Selector class
DS = rmc.Data_Selector(ax, line, proceed_func=stats_funk)

#%%
# ION CALIBRATION A+C
fig_ion = ion_vs_ion(ana, ana.n_sigma)

# HACK to plot multiples figures using the "ion_vs_ion" function
label = 'ION CALIBRATION A+C'
fig_ion.set_label(label)
fig_ion.canvas.set_window_title(label)
fig_ion.set_size_inches(15, 10, forward=True)

axes = fig_ion.get_axes()

for ax in axes:
    ax.set_xlim(-80, 80)
    ax.set_ylim(-80, 80)

# selecting the ionA vs ionC plot
ax = axes[2]
# if lost, modify the title to vizualize the ax in which to make the selection
ax.set_title('')

# HACK get the line corresponding to the quality events
line = ax.lines[1]
assert (line.get_label() == 'Quality Events')

def stats_funk(indexes):
    """
    Used for the proceed button in the Data_selector class.
    Just print the mean and std of the selection along the X and Y axis.
    """
    x_data = line.get_xdata()[indexes]
    y_data = line.get_ydata()[indexes]
    print('X Axis: mean={}    and    std={}'.format(x_data.mean(), x_data.std()))
    print('Y Axis: mean={}    and    std={}'.format(y_data.mean(), y_data.std()))

# calling the Data_Selector class
DS = rmc.Data_Selector(ax, line, proceed_func=stats_funk)

#%%
# actual CALIBRATION of the data (ADU to keV)

# position of the 10.37 keV peak for the different channels
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
        'tk21l000':[371, 1, -55, -55, 55, 55],
        'tk21l001':[304, 1, -55, -55, 55, 55],
        'tk21l002':[266, 1, -55, -55, 55, 55],
        'tk25l000':[567, 1, 50.58, -59, -49.35, 59],
        'tk26l000':[561, 1, 50.58, -59, -51.97, 59],
        'tk26l001':[565, 1, 52, -59, -52.2, 59],
        'tk27l001':[559, 1, 49.15, -59, -49.01, 59],
        'tk27l002':[490, 1, 51, -54, -51, 58],
}


# energy in ADU
energy_raw = ana.all.trig.filt_decor.Energy_OF[ana.all.trig.cut.quality]
energy_raw_nodecor = ana.all.trig.filt.Energy_OF[ana.all.trig.cut.quality]
# energy in keV
energy_calib = energy_raw * 10.37 / np.array(pos10_dict[ana.run])
energy_calib_nodecor = energy_raw_nodecor * 10.37 / np.array(pos10_dict[ana.run])

ion_conv = np.sum(np.array([-1, 1, 1, -1])*energy_calib[:,2:], axis=1)/2
ion_conv_nodecor = np.sum(np.array([-1, 1, 1, -1])*energy_calib_nodecor[:,2:], axis=1)/2


fig_conv, axes = plt.subplots(nrows=2, sharex=True, figsize=(7,8), num='fig_conv')

for ax in axes:
    n, bins, patches = ax.hist(ion_conv_nodecor, bins=250,
            color='deepskyblue', alpha=0.7,
            label='filt\nmean={:.6f}\nstd={:.6f}'.format(ion_conv_nodecor.mean(), ion_conv_nodecor.std())
            )
    ax.hist(ion_conv, bins=bins,
            color='coral', alpha=0.7,
            label='filt decor\nmean={:.6f}\nstd={:.6f}'.format(ion_conv.mean(), ion_conv.std())
            )

    ax.set_ylabel('Quality events counts')
    ax.legend()
    ax.grid()

axes[0].set_title('Charge Conservation')
axes[1].set_xlabel('Total Energy Ionization (A+B+C+D)')
axes[1].set_yscale('log')
fig_conv.tight_layout()
fig_conv.subplots_adjust(hspace=.0)

if save_flag:
    fig_conv.savefig(SAVE_DIR+'/fig_conv.png')

#%%
# SELECTION of the 10keV events
# heat energy
heat = energy_calib[:,0]
# total ionization energy
ion_tot = np.sum(energy_calib[:,2:], axis=1)/2


polar_ion = ana.all.run_tree.Polar_Ion[0]
delta_volt = abs(polar_ion[1]-polar_ion[3])

fig, ax = plt.subplots(num='Tot Ion vs Heat', figsize=(10, 7))
ax.plot(heat, ion_tot, label='quality events', ls='none', marker='.', color='blue', markersize=2)

#guide for 10keV
ax.plot([10.37/(1+delta_volt/3), 10.37], [0, 10.37], 
         zorder=-20, lw=10,
         color='gold', label='10keV band (theory)')

ax.grid()

ax.set_xlim(-2, 13)
ax.set_ylim(-2, 13)

fname_10kev_events = ana.run+'_10kev_events.npy'
fpath_10kev_events = '/'.join((SAVE_DIR, fname_10kev_events))

try:
    crude_ind = np.load(fpath_10kev_events)
    ax.plot(heat[crude_ind], ion_tot[crude_ind],
            label='10.37 keV selection',
            ls='none', marker='.', color='coral', markersize=16,
            zorder=-10)
except:
    print('File not found: {}'.format(fpath_10kev_events))

ax.legend()

def crude_funk(indexes):
    np.save(fpath_10kev_events, indexes)
    print('Indexes were saved in: {}'.format(fpath_10kev_events))

line = ax.lines[0]
DS_2 = rmc.Data_Selector(ax, line, proceed_func=crude_funk)

if save_flag:
    fig.savefig(SAVE_DIR+'/fig_selection.png')


#%%
# ESTIMATING the FIDUCIAL VOLUME
# getting the indexes of the 10keV events
crude_ind = np.load(fpath_10kev_events)
energy_10kev = energy_calib[crude_ind]
time_10kev = ana.all.trig.time[ana.all.trig.cut.quality][crude_ind]

print('Total number: ', crude_ind.size)

fig = plot_ion_vs_ion(ana, energy_10kev, marker='.', markersize=3, alpha=1)

# estimation in the ionB vs ionD plot
axes = fig.get_axes()
ax = axes[3]
ax.set_title('gotcha')

fname_volume = ana.run+'_volume.npy'
fpath_volume = '/'.join((SAVE_DIR, fname_volume))
def count_funk(indexes):
    np.save(fpath_volume, indexes)
    print('Number in Selection: ', indexes.size)

if save_flag:
    fig.savefig(SAVE_DIR+'/ion_10keV.png')

line = ax.lines[0]
DS_3 = rmc.Data_Selector(ax, line, proceed_func=count_funk)

#%%
#heat_adv, _, ionA, ionB, ionC, ionD = energy_10kev.T
ion_list = energy_10kev.T[2:]
cut_list = [ion>2 for ion in ion_list]

time_list = [time_10kev[cut] for cut in cut_list]
ion_cut_list = [ion[cut] for cut,ion in zip(cut_list, ion_list)]

time_mean = list()
ion_med = list()
ion_sup = list()
ion_inf = list()
for ion, time in zip(ion_cut_list, time_list):
    
    npart = int(time[-1]/1.)
    tpart = time[-1] / npart
    
    
    T_mean = list()
    I_med = list()
    I_sup = list()
    I_inf = list()
    for n in range(npart):
        tau = n*tpart
        T_mean.append(tau + tpart/2)
        cut_part = (time > tau) & (time <= tau + tpart)
        ion_part = ion[cut_part]
        
        inf, med, sup = np.quantile(ion_part, [0.16, 0.5, 0.84])
        I_med.append(med)
        I_sup.append(sup-med)
        I_inf.append(med-inf)
        
    time_mean.append(T_mean)
    ion_med.append(I_med)
    ion_inf.append(I_inf)
    ion_sup.append(I_sup)

label_list = ['ion A', 'ion B', 'ion C', 'ion D']

fig, axes = plt.subplots(nrows=4, ncols=2, figsize=(10, 7), num='Histo Ion')

axes_time = axes[:, 0]
axes_hist = axes[:, 1]

for i, ax in enumerate(axes_time):
    
    ion = ion_cut_list[i]
    time = time_list[i]
    label=label_list[i]
    ax.plot(time, ion, ls='none', marker='.', label='10keV events', alpha=0.5)
    
    err = np.array(list((zip(ion_inf[i], ion_sup[i])))).T
    ax.errorbar(time_mean[i], ion_med[i],
               yerr = err,
               ls='none', marker='s', color='r',
               ecolor='k',
               capsize=3,
               capthick=3,
               elinewidth=3,
               label='[0.16, 0.5, 0.84]')
    
    ax.set_ylabel(label)
    ax.legend(loc='lower left')
    ax.grid()

axes_time[-1].set_xlabel('Time [hours]')

for i, ax in enumerate(axes_hist):
    ion = ion_cut_list[i]
    ax.hist(ion ,bins=100, label='10keV events')
    ax.legend(title='median={:.6f}\nstd={:.6f}'.format(np.median(ion), ion.std()),
              loc='upper left')
    ax.grid()

fig.tight_layout()
fig.subplots_adjust(hspace=0., wspace=0.)

if save_flag:
    fig.savefig(SAVE_DIR+'/ion_hist.png')

