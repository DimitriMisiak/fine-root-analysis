#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""
import os
import numpy as np
import matplotlib.pyplot as plt
from spec_classes_fond_neutron import Analysis_red
from representation import (
        temporal_plot, plot_chi2_vs_energy,
        histogram_adu, histogram_ev, ion_vs_ion,
        virtual_vs_virtual_ev, optimization_info,
)
from model_spectrum import fid_mixture, double_norm
from plot_addon import LegendTitle, custom_autoscale, ax_hist, plot_ion_vs_ion, basic_corner
from stats_addon import cdf_calc, custom_bin_edges
import red_magic as rmc


# matplotlib configuration
plt.close('all')
plt.rcParams['text.usetex']=True


stream = 'tg27l000'
save_flag = True

# global variables
SAVE_DIR = '/home/misiak/Analysis/fine_root_analysis/fond_neutron/{}'.format(stream)
DATA_DIR = '/home/misiak/Data/data_run57_neutron/Data'
os.makedirs(SAVE_DIR, exist_ok=True)


thresh_heat = 400
time_cut = None

if stream == 'tg18l005':
    time_cut = (7.4, 7.6)

if stream == 'tg27l000':
    thresh_heat = 700
    time_cut = (7, 11.3)

if stream == 'tg28l000':
    thresh_heat = 700
    time_cut = (7.4, 8.05)

thresh_ion = 300
# getting the data
ana = Analysis_red(
        stream,
        detector='RED80',
        run_dir=DATA_DIR,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=5,
        override_mu=(861, -54, 54),
        override_sigma=(40, 10, 9),
        thresh_heat=thresh_heat,
        thresh_ion=thresh_ion,
        time_cut= time_cut,
)

# =============================================================================
# PLOTS
# =============================================================================

def line_in_new_ax(ax, line):
    x_data, y_data = line.get_xdata(), line.get_ydata()
    color = line.get_color()
    ls = line.get_linestyle()
    lw = line.get_linewidth()
    zorder = line.get_zorder()
    label = line.get_label()
    marker = line.get_marker()
    alpha =line.get_alpha()
    
    new_lines = ax.plot(
            x_data, y_data,
            color=color, ls=ls, lw=lw, marker=marker, alpha=alpha,
            zorder=zorder, label=label, 
    )
    
    return new_lines

# monitoring vs TIME
fig_temp = temporal_plot(ana)

# heat vsv Time
fig_heat , ax = plt.subplots(figsize=(10, 7), num='Heat vs Time')
for line in fig_temp.get_axes()[0].lines:
    line_in_new_ax(ax, line)
ax.set_ylabel('Energy heat [ADU]')
ax.set_xlabel('Time [hours]')
ax.set_yscale('log')
ax.set_ylim(10**1, 10**5)


# QUALITY CUTS
fig_chi2_trig, fig_chi2_noise = plot_chi2_vs_energy(ana)

# plotting the quality cuts
x_data = 10**np.linspace(-2, 5, int(1e4))
cut_ion = thresh_ion*(1+(x_data/3e2)**2.2)
cut_heat = thresh_heat*(1+(x_data/2e3)**2)
for i, ax in enumerate(fig_chi2_trig.get_axes()):
    if i == 2:
        ax.plot(x_data, cut_heat, lw=1, color='k', label='quality cut')
    else:
        ax.plot(x_data, cut_ion, lw=1, color='k', label='quality cut')        
        

# resize the plots
for ax in fig_chi2_trig.get_axes():
    ax.set_xlim(10**-2, 10**5)
    ax.set_ylim(10**1, 10**9)
    ax.legend()

# HEAT CALIBRATION
fig_hist_trig, fig_hist_noise = histogram_adu(ana)

# resize the plots
fig_hist_trig.get_axes()[0].set_xlim(-200, 2000)
for i, ax in enumerate(fig_hist_trig.get_axes()[:5]):
    if i==0:
        ax.set_xlim(-200, 2000)
    else:
        print(i)
        ax.set_xlim(-70, 70)
    
    
    
if save_flag:
    fig_temp.savefig(SAVE_DIR+'/fig_temp.png')
    fig_heat.savefig(SAVE_DIR+'/fig_heat.png')
    fig_chi2_trig.savefig(SAVE_DIR+'/fig_chi2_trig.png')
    fig_chi2_noise.savefig(SAVE_DIR+'/fig_chi2_noise.png')
    fig_hist_trig.savefig(SAVE_DIR+'/fig_hist_trig.png')
    fig_hist_noise.savefig(SAVE_DIR+'/fig_hist_noise.png')   

    
#%%
# CROSSTALK correction
energy_array = ana.all.trig.filt_decor.Energy_OF[ana.all.trig.cut.quality]
#energy_array = ana.all.trig.filt.Energy_OF[ana.all.trig.cut.quality]

#corr_matrix = np.array([
#        [1, ab, ac, ad],
#        [ba, 1, bc, bd],
#        [ca, cb, 1, cd],
#        [da, db, dc, 1]
#])
fig_cross = plot_ion_vs_ion(ana, energy_array, alpha=0.1, zorder=-5)
axes = fig_cross.get_axes()
for ax in axes:
    ax.axvline(0, color='r')
    ax.axhline(0, color='r')

A, B, C, D = energy_array[:, 2:].T
line0, = axes[0].plot(D, A, marker='.', ls='none', markersize=3)
line1, = axes[1].plot(D, C, marker='.', ls='none', markersize=3)
line2, = axes[2].plot(A, C, marker='.', ls='none', markersize=3)
line3, = axes[3].plot(D, B, marker='.', ls='none', markersize=3)
line4, = axes[4].plot(A, B, marker='.', ls='none', markersize=3)
line5, = axes[5].plot(C, B, marker='.', ls='none', markersize=3)



#%%
#corr_matrix = np.array([
#        [1, 0, 0, 0.05],
#        [0, 1, 0, 0],
#        [0, 0, 1, -0.025],
#        [0.05, -0.02, 0, 1]
#])

corr_matrix = np.array([
        [1, -0.052, 0, 0],
        [-0.03, 1, 0, 0],
        [-0.012, 0.001, 1, -0.025],
        [0, 0, -0.03, 1]
])
    
A, B, C, D = np.dot(corr_matrix, energy_array[:, 2:].T)

line0.set_xdata(D)
line0.set_ydata(A)

line1.set_xdata(D)
line1.set_ydata(C)

line2.set_xdata(A) 
line2.set_ydata(C)
   
line3.set_xdata(D)
line3.set_ydata(B)

line4.set_xdata(A)
line4.set_ydata(B)

line5.set_xdata(C)
line5.set_ydata(B)

if save_flag:
    fig_cross.savefig(SAVE_DIR+'/fig_cross.png')

#%%

ana.all.trig.filt_decor.Energy_OF[:, 2:] = np.dot(corr_matrix, ana.all.trig.filt_decor.Energy_OF[:, 2:].T).T

ana.all.trig.filt.Energy_OF[:, 2:] = np.dot(corr_matrix, ana.all.trig.filt.Energy_OF[:, 2:].T).T


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
ax = axes[2]
#ax = axes[3]
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
        'tg17l007':[1165, 1, -56, -54, 57, 56],
        'tg18l005':[1177, 1, -56, -54, 57, 56],
        'tg19l010':[1199, 1, -56, -54, 57, 56],
        'tg20l000':[1197, 1, -56, -54, 57, 56],
        'tg21l000':[1235, 1, -56, -54, 57, 56],
        'tg27l000':[1293, 1, -56, -54, 57, 56],
        'tg28l000':[1291, 1, -56, -54, 57, 56],
}


# energy in ADU
energy_raw = ana.all.trig.filt_decor.Energy_OF[ana.all.trig.cut.quality]
energy_raw_nodecor = ana.all.trig.filt.Energy_OF[ana.all.trig.cut.quality]
# energy in keV


energy_calib = energy_raw * 10.37 / np.array(pos10_dict[ana.run])
energy_calib_nodecor = energy_raw_nodecor * 10.37 / np.array(pos10_dict[ana.run])
### HACK
#energy_calib = energy_calib_nodecor

ion_conv = np.sum(np.array([1, 1, -1, -1])*energy_calib[:,2:], axis=1)/2
ion_conv_nodecor = np.sum(np.array([1, 1, -1, -1])*energy_calib_nodecor[:,2:], axis=1)/2

    
threshold_conv = 1
cut_conv = (abs(ion_conv)<=threshold_conv)
cut_noconv = (abs(ion_conv)>threshold_conv)

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


    ax.axvline(threshold_conv, color='k', ls='--')
    ax.axvline(-threshold_conv, color='k', ls='--')
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
# charge conservation vs each ionization channel
    
#samples = np.c_[energy_calib_nodecor[:, 2:], ion_conv]
samples = np.c_[energy_calib[:, 2:], ion_conv]

labels = ('ionA', 'ionB', 'ionC', 'ionD', 'Conv')
fig, axes = basic_corner(samples[cut_conv, :], labels, color='k')
fig, axes = basic_corner(samples[cut_noconv, :], labels, axes=axes, color='r')

#%%
# SELECTION of the 10keV events
# heat energy
heat = energy_calib[:,0]
# total ionization energy
ion_tot = np.sum(energy_calib[:,2:], axis=1)/2
#ion_tot = energy_calib[:, 5]


polar_ion = ana.all.run_tree.Polar_Ion[0]
delta_volt = abs(polar_ion[1]-polar_ion[3])

fig, ax = plt.subplots(num='Tot Ion vs Heat', figsize=(10, 7))
ax.plot(heat[cut_conv], ion_tot[cut_conv], label='quality events', ls='none', marker='.', color='k', markersize=2)
ax.plot(heat[cut_noconv], ion_tot[cut_noconv], label='quality events', ls='none', marker='.', color='r', markersize=2)

#guide for 10keV
ax.plot([10.37/(1+delta_volt/3), 10.37], [0, 10.37], 
         zorder=-20, lw=10,
         color='gold', label='10keV band (theory)')

ax.grid()

ax.set_xlim(-2, 13)
ax.set_ylim(-2, 13)
ax.set_ylabel('Total Ionization Energy A+B+C+D [keV]')
ax.set_xlabel('Heat Energy [keV]')
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
#DS_2 = rmc.Data_Selector(ax, line, proceed_func=crude_funk)

if save_flag:
    fig.savefig(SAVE_DIR+'/fig_selection.png')

