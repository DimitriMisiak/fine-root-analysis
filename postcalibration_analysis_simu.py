#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 09:17:47 2019

@author: misiak
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import json
from plot_addon import basic_corner
import red_magic as rmc
import matplotlib.patheffects as pe

# matplotlib configuration
plt.close('all')
plt.rcParams['text.usetex']=True
cartoon = [
        pe.Stroke(linewidth=3, foreground='k'),
        pe.Normal(),
]


simu_list = [
        'Flat_Analytical_SimuOnly_0.0000_50.0000_ER',
        'Flat_Analytical_SimuOnly_0.0000_50.0000_NR',
        'Line_Analytical_SimuOnly_1.3000_ER',
        'Line_Analytical_SimuOnly_10.3700_ER',
]

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


#stream = 'tg27l000'
simu = simu_list[2]
save_flag = True

with open('/home/misiak/Data/data_run57_neutron/stream_config.json', 'r') as cfg:
    config = json.load(cfg)

calibration_quality = list()
for stream in config['Calibration']:
    save_dir = '/home/misiak/Analysis/fine_root_analysis/fond_neutron/{}/{}_tot'.format(stream, simu)
    raw_load = np.load(save_dir+'/energy_quality.npy')
    data_quality = np.delete(raw_load, 1, axis=1)
    calibration_quality.append(data_quality)

energy_quality = np.concatenate(calibration_quality)

#energy_all = np.load(SAVE_DIR+'/energy_calibrated.npy')
#energy_quality = np.load(SAVE_DIR+'/energy_quality.npy')
#
## removing the heatB channel
#energy_all = np.delete(energy_all, 1, axis=1)
#energy_quality = np.delete(energy_quality, 1, axis=1)


heat, A, B, C, D = energy_quality.T

# =============================================================================
# FIDUCIAL CUT
# =============================================================================

fid_thresh = 0.5
# HACK THE POSITIVE HEAT IN and the limit to 50keV
fid_cut = (abs(A)<fid_thresh) & (abs(C)<fid_thresh) & (heat > 0.025) & (heat < 50)


fig, ax = plt.subplots()

ax.plot(A, C,
        label='quality events',
        ls='none', marker='.', color='k', alpha=0.3)

ax.plot(C[fid_cut], A[fid_cut],
        label='fiducial events',
        ls='none', marker='.', color='r', alpha=0.3)

ax.axhspan(-fid_thresh, fid_thresh,
           alpha=0.3, color='crimson', zorder=-5,
           label='fiducial cut')
ax.axvspan(-fid_thresh, fid_thresh,
           alpha=0.3, color='crimson', zorder=-5)

ax.set_xlim(-2, 12)
ax.set_ylim(-2, 12)

ax.set_ylabel('Ion C')
ax.set_xlabel('Ion A')
ax.grid()
ax.legend()


labels = ('heatA', 'ionA', 'ionB', 'ionC', 'ionD')
fig, axes = basic_corner(energy_quality, labels, color='k')
fig, axes = basic_corner(energy_quality[fid_cut], labels, axes=axes, color='r')

#%%
# =============================================================================
# # NR CUT
# =============================================================================
def energy_recoil(ec, ei, V):
#    coeff = 1.6e-19 * V / 3
    coeff = V / 3
    return ec*(1+coeff) - ei*coeff

def quenching(ec, ei, V):
    er = energy_recoil(ec, ei, V)
    return ei/er

def lindhard(er):
    
    A = 72.63
    Z = 32
    
    k = 0.133 * Z**(2./3) * A**(-1./2)
    epsilon = 11.5 * er * Z**(-7./3)
    g = 3 * epsilon**0.15 + 0.7 * epsilon**0.6 + epsilon  
    
    Q = k*g/(1+k*g)
    
    return Q

heat_fid, A_fid, B_fid, C_fid, D_fid = energy_quality[fid_cut].T

collect_fid = (B_fid + D_fid)/2
collect = (B + D)/2

fig, ax = plt.subplots()

#ax.plot(heat, collect,
#        label='quality events',
#        ls='none', marker='.', alpha=0.3, color='grey')

line, = ax.plot(heat_fid, collect_fid,
        label='fid events',
        ls='none', marker='.', alpha=0.1, color='grey')

## calling the Data_Selector class
#DS = rmc.Data_Selector(ax, line, proceed_func=stats_funk)

er_array = np.linspace(0, 200, int(1e4))
dv=2
ec_array = er_array * (1 + lindhard(er_array)*dv/3) / (1 + dv/3)
ei_array = er_array * lindhard(er_array)

#ax.plot(ec_array, ei_array, label='lindhard')

#lindhard_from_ec = np.interp()

#ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_xlim(-0.2, 100)
ax.set_xlabel('Heat Energy[keV]')
ax.set_ylabel('Collect Ionization Energy (B+D)/2 [keV]')
ax.grid()
ax.legend()



#nblob_cut = (heat_fid > 0.05) & (heat_fid < 0.4)
#ax.plot(heat_fid[nblob_cut], collect_fid[nblob_cut],
#        label='noise',
#        ls='none', marker='.', alpha=1, zorder=15, color='k')
#std_nb = np.std(collect_fid[nblob_cut])
std_nb = 0.254
#std_1kev = 0.272
std_10kev = 0.318

alpha = (std_10kev**2 - std_nb**2)**0.5 / 10.37
def std_collect(ec):
    return ( std_nb**2 + (alpha*ec)**2 )**0.5

# gamma cut
ec_graph = np.linspace(0, 100, int(1e4))
std_graph = std_collect(ec_graph)
ax.plot(ec_graph, ec_graph, 
        color='deepskyblue',  path_effects=cartoon,
        zorder=20)
ax.fill_between(ec_graph, ec_graph + 3*std_graph, ec_graph - 3*std_graph,
                color='deepskyblue', alpha=0.2,
                zorder=-10)

std_fid = std_collect(heat_fid)
gamma_cut = (collect_fid < (heat_fid + 3*std_fid)) & (collect_fid > (heat_fid - 3*std_fid))
ax.plot(heat_fid[gamma_cut], collect_fid[gamma_cut],
        label='gamma band',
        ls='none', marker='.', alpha=0.7, markersize=10, color='steelblue')

# neutron cut
ei_graph = np.interp(ec_graph, ec_array, ei_array)
ax.plot(ec_graph, ei_graph, 
        color='coral',  path_effects=cartoon,
        zorder=20)
ax.fill_between(ec_graph,  ei_graph + 3*std_graph, ei_graph - 3*std_graph,
                color='coral', alpha=0.2,
                zorder=-10)

collect_lindhard = np.interp(heat_fid, ec_array, ei_array)
neutron_cut = (collect_fid < (collect_lindhard + 3*std_fid)) & (collect_fid > (collect_lindhard - 3*std_fid))
ax.plot(heat_fid[neutron_cut], collect_fid[neutron_cut],
        label='neutron band',
        ls='none', marker='.', alpha=0.7, markersize=10, color='coral')

ax.legend()

# HO cut
ax.plot(ec_graph, np.zeros(ec_graph.shape), 
        color='k',  path_effects=cartoon,
        zorder=20)
ax.fill_between(ec_graph,  3*std_graph, - 3*std_graph,
                color='k', alpha=0.2,
                zorder=-10)

HO_cut = (collect_fid < (0 + 3*std_fid)) & (collect_fid > (0 - 3*std_fid))
ax.plot(heat_fid[HO_cut], collect_fid[HO_cut],
        label='HO band',
        ls='none', marker='.', alpha=0.7, markersize=10, color='k', zorder=-10)

ax.legend()



#neutron_cut = np.abs(collect_fid - )

### 1 KEV SECTION
#heat_start = 1
#heat_window = 1
#
#heat_bar = list()
#collect_bar = list()
#for i in range(30):
#    
#    heat_inf = heat_start + heat_window*i
#    heat_sup = heat_inf + heat_window
#    section_cut = (heat_fid > heat_inf) & (heat_fid < heat_sup)
#    
#    heat_section = heat_fid[section_cut]
#    collect_section = collect_fid[section_cut]
#    
#    heat_bar.append( np.quantile(heat_section, [0.16, 0.5, 0.84]) )
#    collect_bar.append( np.quantile(collect_section, [0.16, 0.5, 0.84]) )
#    
#heat_bar_array = np.array(heat_bar)
#collect_bar_array = np.array(collect_bar)
#
#heat_inf, heat_med, heat_sup = heat_bar_array.T
#collect_inf, collect_med, collect_sup = collect_bar_array.T
#sig_inf = collect_med - collect_inf
#sig_sup = collect_sup - collect_med
#
#ax.errorbar(heat_med, collect_med, yerr=(sig_inf, sig_sup),
#            ls='none', marker='o', color='k')


#%%

# =============================================================================
# QUENCHING PLOT
# =============================================================================
dv = 2 #2V
er = energy_recoil(heat_fid, collect_fid, dv)
q = quenching(heat_fid, collect_fid, dv)

er_array = np.linspace(0, 50, int(1e4))
quenching_array = 0.16*(er_array)**0.18


fig, ax = plt.subplots()

#ax.plot(heat, collect,
#        label='quality events',
#        ls='none', marker='.', alpha=0.3, color='grey')

ax.plot(er, q,
        label='fid events',
        ls='none', marker='.', color='grey', alpha=0.3, zorder=-10)

ax.plot(er[gamma_cut], q[gamma_cut],
        label='in gamma band',
        ls='none', marker='.', color='steelblue', alpha=0.3)

ax.plot(er[neutron_cut], q[neutron_cut],
        label='in NR band',
        ls='none', marker='.', color='coral', alpha=0.3)

ax.plot(er[HO_cut], q[HO_cut],
        label='in HO band',
        ls='none', marker='.', color='k', alpha=0.3, zorder=-5)


#ax.plot(er_array, quenching_array,
#        label='NR quenching theory')

ax.plot(er_array, lindhard(er_array),
        label='Lindhard')

ax.set_xlim(0, 30)
ax.set_ylim(-0.25, 1.25)
ax.set_xlabel('Recoil Energy [keV]')
ax.set_ylabel('Quenching factor')
ax.grid()
ax.legend()

#%%
# =============================================================================
# Energy Spectrum
# =============================================================================

# keVee
fig, axes = plt.subplots(ncols=2, figsize=(10, 5))
fig.suptitle('Heat Energy Spectrum keVee')

for ax in axes:
    n, bins, patches = ax.hist(heat_fid, bins=250,
                               label='fid events',
                               color='grey', alpha=0.1)
    ax.hist(heat_fid[gamma_cut], bins=bins,
            label='Gamma band',
            color='deepskyblue', alpha=0.5)
    ax.hist(heat_fid[neutron_cut], bins=bins,
            label='NR band',
            color='coral', alpha=0.5)
#    ax.hist(heat_fid[HO_cut], bins=bins,
#            label='HO band',
#            color='k', alpha=0.1)
    
    ax.set_ylabel('Counts, bin width = {:.3f}'.format(bins[1]-bins[0]))
    ax.set_xlabel('Heat Energy [keVee]')
    ax.set_yscale('log')
    ax.legend()
axes[1].set_xscale('log')


# keV
heat_gamma = energy_recoil(heat_fid[gamma_cut], collect_fid[gamma_cut], dv)
heat_neutron = energy_recoil(heat_fid[neutron_cut], collect_fid[neutron_cut], dv)
heat_ho = energy_recoil(heat_fid[HO_cut], collect_fid[HO_cut], dv)

fig, axes = plt.subplots(ncols=2, figsize=(10, 5))
fig.suptitle('Recoil Energy Spectrum keV')

for ax in axes:

#n, bins, patches = ax.hist(heat_fid, bins=250,
#                           label='fid events',
#                           color='grey', alpha=0.1)
    
    n, bins, patches = ax.hist(
            heat_gamma, bins=500,
            label='Gamma band',
            color='deepskyblue', alpha=0.5
    )
    ax.hist(heat_neutron, bins=bins,
            label='NR band',
            color='coral', alpha=0.5)
#    ax.hist(heat_ho, bins=bins,
#            label='HO band',
#            color='k', alpha=0.1)

    ax.set_ylabel('Counts, bin width = {:.3f}'.format(bins[1]-bins[0]))
    ax.set_xlabel('Recoil Energy [keV]')
    ax.set_yscale('log')
    ax.legend()
axes[1].set_xscale('log')
