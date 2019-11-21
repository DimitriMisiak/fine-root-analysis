#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""

import numpy as np
import matplotlib.pyplot as plt
from spec_classes import Analysis_red
from representation import (
        temporal_plot, plot_chi2_vs_energy,
        histogram_adu, histogram_ev, ion_vs_ion,
        virtual_vs_virtual_ev, optimization_info
)
from model_spectrum import fid_mixture, double_norm
from plot_addon import LegendTitle, custom_autoscale, ax_hist
from stats_addon import cdf_calc, custom_bin_edges

run_dir = '/home/misiak/Data/data_run57'
              

detector = 'RED80'


# first command
plt.close('all')
plt.rcParams['text.usetex']=True

n_sigma = 2

### for tg12l003
stream = 'tg12l003'  
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

### for tg13l000
stream = 'tg13l000' 
ana_2 = Analysis_red(
        stream,
        detector=detector,
        run_dir=run_dir,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=n_sigma,
        override_mu=(900, -53.4, 55),
        override_sigma=(60, 10, 8)
)

### for tg14l001
stream = 'tg14l001'
ana_3 = Analysis_red(
        stream,
        detector=detector,
        run_dir=run_dir,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=n_sigma,
        override_mu=(928, -53, 54),
        override_sigma=(40, 8, 9)
)   

### for tg19l010
stream = 'tg19l010' 
ana_4 = Analysis_red(
        stream,
        detector=detector,
        run_dir=run_dir,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=n_sigma,
        override_mu=(1200, -54.4, 55.7),
        override_sigma=(50, 7, 8)
)

### for tg20l000
stream = 'tg20l000' 
ana_5 = Analysis_red(
        stream,
        detector=detector,
        run_dir=run_dir,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=n_sigma,
        override_mu=(1200, -54.4, 55.7),
        override_sigma=(60, 7, 8)
)

### for tg21l000
stream = 'tg21l000'
ana_6 = Analysis_red(
        stream,
        detector=detector,
        run_dir=run_dir,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0, 3, 5),
        n_sigma=n_sigma,
        override_mu=(1236, -54.4, 55.7),
        override_sigma=(50, 7, 8)
)   

#%%

plt.close('all')

##### USUAL PLOTSSS
#ana=ana_6
#fig_temp = temporal_plot(ana)
#fig_chi2_trig, fig_chi2_noise = plot_chi2_vs_energy(ana)
#fig_hist_trig, fig_hist_noise = histogram_adu(ana)
#fig_hist_trig_ev, fig_hist_noise_ev = histogram_ev(ana)
#fig_ion = ion_vs_ion(ana, n_sigma)
#fig_virtual = virtual_vs_virtual_ev(ana)

#### ALL STREAMS PLOTS
ana_list = [ana_1, ana_2, ana_3, ana_4, ana_5, ana_6]

stream_list = [ana.run for ana in ana_list]
sigma0_heat_list = [ana.all.noise.sigma0_ev.heat_a for ana in ana_list]
sigma0_collect_list = [ana.all.noise.sigma0_ev.collect for ana in ana_list]
sigma0_collect_list[1] = 246

energy_list = [ana.all.trig.energy_ev for ana in ana_list]

cut_fid_list = [ana.all.trig.cut.fiducial for ana in ana_list]
cut_qual_list = [ana.all.trig.cut.quality for ana in ana_list]

heat_qual_list = [energy.heat_a[cut] for energy, cut in zip(energy_list, cut_qual_list)]
heat_fid_list = [energy.heat_a[cut] for energy, cut in zip(energy_list, cut_fid_list)]
collect_qual_list = [energy.collect[cut] for energy, cut in zip(energy_list, cut_qual_list)]
collect_fid_list = [energy.collect[cut] for energy, cut in zip(energy_list, cut_fid_list)]

heat_qual = np.concatenate(heat_qual_list)
heat_fid = np.concatenate(heat_fid_list)
collect_qual = np.concatenate(collect_qual_list)
collect_fid = np.concatenate(collect_fid_list)

def energy_recoil(ec, ei, V):
#    coeff = 1.6e-19 * V / 3
    coeff = V / 3
    return ec*(1+coeff) - ei*coeff

def quenching(ec, ei, V):
    er = energy_recoil(ec, ei, V)
    return ei/er

dv = 2 #2V
er_qual = energy_recoil(heat_qual, collect_qual, dv)
quenching_qual = quenching(heat_qual, collect_qual, dv)
er_fid = energy_recoil(heat_fid, collect_fid, dv)
quenching_fid = quenching(heat_fid, collect_fid, dv)

#bin_size = 1000
#bins_energy_inf = np.arange(0, 12000, bin_size)
#bins_energy_sup = np.arange(0, 12000, bin_size) + bin_size

#bins_energy_inf = [0, 1000, 3000, 4000, 6000, 7500, 10000]
#bins_energy_sup = [1000, 2000, 3500, 5000, 6500, 8000, 11000]

bins_energy_inf = [700, 1000, 1300, 3000, 10000]
bins_energy_sup = [900,1200, 1500, 3500, 11000]

d_norm = double_norm()
heat_binned = list()
collect_binned = list()

for inf, sup in zip(bins_energy_inf, bins_energy_sup):
    truth_inf = heat_fid > inf
    truth_sup = heat_fid < sup
    truth_tot = np.logical_and(truth_inf, truth_sup)
    heat_binned.append(heat_fid[truth_tot])
    collect_binned.append(collect_fid[truth_tot])
    

popt_list = [d_norm.fit(collect, floc=0, fscale=1) for collect in collect_binned]

# histo resolution
fig = plt.figure()
plt.bar(stream_list, sigma0_heat_list, width=0.4, align='edge', label='Heat')
plt.bar(stream_list, sigma0_collect_list, width=0.3, align='center', label='Ion B+D')
plt.xlabel('Stream name')
plt.ylabel('Resolution [eV]')
plt.legend()
plt.grid()


fig, axes = plt.subplots(nrows=2, figsize=(10,10))

for ax in axes:
    
    ax.plot(
            heat_fid, collect_fid,
            ls='none', marker='2', zorder=11, color='slateblue',
            label='Fiducial Events', alpha=0.7
    )
    ax.plot(
            heat_qual, collect_qual,
            ls='none', marker='1', zorder=10, color='coral',
            label='Quality Events', alpha=0.7
    )
      
    ax.legend(title='All streams')
    ax.set_xlabel('Heat Energy [eV]')
    ax.set_ylabel('Ion Energy [eV]')
    ax.grid(alpha=0.3)

    
axes[0].set_xlim(0, 12e3)
axes[1].set_xlim(8e2, 12e3)
axes[0].set_ylim(-1000, 12e3)
axes[1].set_ylim(100, 20e3)
axes[1].set_xscale('log')
axes[1].set_yscale('log')

fig.tight_layout()

####### for quenching vs ER

er_array = np.linspace(0, 15e3, int(1e4))
quenching_array = 0.16*(er_array*1e-3)**0.18

fig, axes = plt.subplots(nrows=1, figsize=(10,10))
axes = [axes,]
for ax in axes:
    
    ax.plot(
            er_qual, quenching_qual,
            ls='none', marker='1', zorder=10, color='coral',
            label='Quality Events', alpha=0.7
    )    
    
    ax.plot(
            er_array, quenching_array,
            lw=5, zorder=10, color='limegreen',
            label='$Q=0.16E_R^{0.18}$', alpha=0.7
    )    
        
    ax.plot(
            er_fid, quenching_fid,
            ls='none', marker='2', zorder=11, color='slateblue',
            label='Fiducial Events', alpha=1,
    )

      
    ax.legend(title='All streams')
    ax.set_xlabel('Recoil Energy [eV]')
    ax.set_ylabel('Quencing factor')
    ax.grid(alpha=0.3)

fig.tight_layout()








####### for section histogrmm

#fig, axes = plt.subplots(nrows=1, figsize=(10,10))
#axes = [axes,]
#for ax in axes:
#    
#    ax.plot(
#            heat_fid, collect_fid,
#            ls='none', marker='o', zorder=11, color='slateblue',
#            label='Fiducial Events', alpha=0.1
#    )
##    ax.plot(
##            heat_qual, collect_qual,
##            ls='none', marker='1', zorder=10, color='coral',
##            label='Quality Events', alpha=0.7
##    )
#      
#    ax.legend(title='All streams')
#    ax.set_xlabel('Heat Energy [eV]')
#    ax.set_ylabel('Ion Energy [eV]')
#    ax.grid(alpha=0.3)
#
#    for inf, sup in zip(bins_energy_inf, bins_energy_sup):
#        ax.axvspan(inf, sup, alpha=0.7, color='coral')
#    
#axes[0].set_xlim(0, 12e3)
##axes[1].set_xlim(8e2, 12e3)
#axes[0].set_ylim(-1000, 12e3)
##axes[1].set_ylim(100, 20e3)
##axes[1].set_xscale('log')
##axes[1].set_yscale('log')
#
#fig.tight_layout()


#
## section hist
#nbin = len(heat_binned)
#
##fig, axes = plt.subplots(ncols=int(nbin**0.5)+1, nrows=int(nbin**0.5), figsize=(15, 10))
##axes = axes.flatten()
#for i in range(nbin):
#    
#    fig, ax = plt.subplots()
#    
#    
#    heat = heat_binned[i]
#    collect = collect_binned[i]
#    popt = popt_list[i]
#    if collect.size == 0:
#        continue
#    
##    ax = axes[i]
#
##    ax.hist(collect, bins=50)
#    bin_edges = custom_bin_edges(collect, 50)
#        
#    a0 = ax_hist(ax, bin_edges, collect,
#            'Fiducial events', color='limegreen')[0]   
#    
#
#    xrange = np.linspace(collect.min(), collect.max(), 1000)
#    pdf = d_norm.pdf(xrange, *popt)
#    cdf = d_norm.cdf(xrange, *popt)
#    normalization = collect.size
#    pdf_norm = pdf * normalization * (bin_edges[1] - bin_edges[0])
#    
#    ax.autoscale(False)
#    ax.plot(xrange, pdf_norm,
#            ls='--', color='yellow',
#            label='fit')
#    
#    a0.plot(xrange, cdf,
#            ls='-.', color='yellow',
#            label='fit')
#