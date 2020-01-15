#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:04:43 2020

@author: misiak
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from representation_functions import (
    temporal_plot,
    temporal_plot_heat_only,
    plot_chi2_vs_energy,
    histogram_adu,
    histogram_ev,
    ion_vs_ion,
    virtual_vs_virtual_ev
)

from plot_addon import basic_corner

from data_analysis import (
    analysis_parameters,
    heat_chi2_threshold_function,
    ion_chi2_threshold_function,
)


plt.close('all')
plt.rcParams['text.usetex']=True

analysis_dir = '/home/misiak/Analysis/neutron_background'
analysis_data_path = '/'.join([analysis_dir, 'simu_analysis.h5'])

aux_df = pd.read_hdf(
    analysis_data_path,
    key='df',
    columns=['stream', 'simulation']
)
stream_list = aux_df['stream'].unique()    
simulation_list = aux_df['simulation'].unique()   

config_list = list()
for stream in stream_list:
    for simu in simulation_list:
        config_list.append( [stream, simu] )

from tqdm import tqdm
for config in tqdm(config_list):
# if True:
#     config = ['tg28l000', 'line_10keV']    
    
    stream, simulation = config
    
    plt.close('all')
    
    df_analysis = pd.read_hdf(
        analysis_data_path,
        key='df',
        where='stream = "{}" & simulation = "{}"'.format(stream, simulation)
    )
    
    source = df_analysis['source'].unique()[0]    
    
    ### monitoring
    fig_temp = temporal_plot(stream, df_analysis)
    fig_temp_heat = temporal_plot_heat_only(stream, df_analysis)
    
    ### chi2 plot
    fig_chi2 = plot_chi2_vs_energy(stream, df_analysis)
    # plotting the quality cuts
    x_data = 10**np.linspace(-2, 5, int(1e4))
    cut_ion = ion_chi2_threshold_function(
        analysis_parameters['ion_chi2_threshold'],
        x_data
    )
    cut_heat = heat_chi2_threshold_function(
        analysis_parameters[stream]['heat_chi2_threshold'],
        x_data)    
    for i, ax in enumerate(fig_chi2.get_axes()):
        if i == 2:
            ax.plot(x_data, cut_heat, lw=1, color='k', label='quality cut')
        else:
            ax.plot(x_data, cut_ion, lw=1, color='k', label='quality cut')
        ax.set_xlim(10**-2, 10**5)
        ax.set_ylim(10**1, 10**9)
        ax.legend()            

    ### histogramm ADU
    fig_hist_trig = histogram_adu(stream, df_analysis, bins=10000)
    # resize the plots
    fig_hist_trig.get_axes()[0].set_xlim(-200, 2000)
    for i, ax in enumerate(fig_hist_trig.get_axes()[:5]):
        if i==0:
            ax.set_xlim(-200, 2000)
        else:
            ax.set_xlim(-70, 70)    

    ### crosstalk correction
    samples = df_analysis[df_analysis.quality_cut][[
        'energy_adu_ionA',
        'energy_adu_ionB',
        'energy_adu_ionC',
        'energy_adu_ionD',
    ]]
    samples_corr = df_analysis[df_analysis.quality_cut][[
        'energy_adu_corr_ionA',
        'energy_adu_corr_ionB',
        'energy_adu_corr_ionC',
        'energy_adu_corr_ionD',
    ]]
    fig_cross, axes = basic_corner(
        samples.values,
        samples.columns,
        num = '{} ({}): Cross-talk Correction'.format(stream, source),
        label='raw',
        alpha=0.1,
    )
    basic_corner(
        samples_corr.values,
        samples_corr.columns,
        axes=axes,
        color='slateblue',
        zorder=-1,
        label='corrected'
    )
    for ax in fig_cross.get_axes():
        ax.axvline(0, color='r', zorder=-5)
        ax.axhline(0, color='r', zorder=-5)

    ### histogramm kev
    fig_hist_trig_ev = histogram_ev(stream, df_analysis, bins=10000)
    for ax in fig_hist_trig_ev.get_axes():
        ax.set_xlim(-2.5, 15)
    
    ### ion vs ion
    fig_ion = ion_vs_ion(stream, df_analysis)
    axes = fig_ion.get_axes()
    for ax in axes:
        ax.set_xlim(-15, 15)
        ax.set_ylim(-15, 15)
    
    ### ion vs ion virtual
    fig_virtual = virtual_vs_virtual_ev(stream, df_analysis)
    for ax in fig_virtual.get_axes():
        ax.set_xlim(-15, 15)
        ax.set_ylim(-15, 15)
    
    ### 10kev plot
    delta_volt = 2 #V
    quality_cut = df_analysis['quality_cut']
    fig_10kev, ax = plt.subplots(num='Tot Ion vs Heat', figsize=(10, 7))
    ax.set_title('{} ({}): 10keV events'.format(stream, source))
    ax.plot(
        df_analysis[quality_cut]['energy_heat'],
        df_analysis[quality_cut]['energy_ion_total'],
        label='quality events',
        ls='none',
        marker='.',
        color='b',
        alpha=0.3
    )
    #guide for 10keV
    ax.plot(
        [10.37/(1+delta_volt/3), 10.37],
        [0, 10.37], 
        zorder=-20, lw=10,
        color='gold', label='10keV band (theory)'
    )
    ax.grid()
    ax.set_xlim(-2, 13)
    ax.set_ylim(-2, 13)
    ax.set_ylabel('Total Ionization Energy A+B+C+D [keV]')
    ax.set_xlabel('Heat Energy [keV]')
    fig_10kev.tight_layout()

    # trigger cut (specific to simulation)
    trigger_cut = df_analysis['trigger_cut']
    fig_trigger, axes = plt.subplots(nrows=2, sharex=True)
    axes[0].plot(
        df_analysis['timestamp'],
        df_analysis['t_input_simu_trigger'],
        label='all simu',
        ls='none', marker='.', markersize=3
    )
    axes[0].plot(
        df_analysis[trigger_cut]['timestamp'],
        df_analysis[trigger_cut]['t_input_simu_trigger'],
        label='trig simu',
        ls='none', marker='.', color='k'
    )
    axes[0].axhspan(-5, 5, color='limegreen')
    axes[1].plot(
            df_analysis['timestamp'],
            df_analysis['t_nearest_data_trigger'],
            label='all simu',
            ls='none', marker='.', markersize=3
    )
    axes[1].plot(
            df_analysis[trigger_cut]['timestamp'],
            df_analysis[trigger_cut]['t_nearest_data_trigger'],
            label='trig simu',
            ls='none', marker='.', color='k'
    )
    axes[1].axhspan(-5, 5, color='coral')
    axes[0].set_ylabel(r'$\Delta T_{simu}$')
    axes[1].set_ylabel(r'$\Delta T_{data}$')
    axes[1].set_xlabel('Time [hours]')
    axes[0].legend()
    axes[1].legend()
    fig_trigger.tight_layout()
    fig_trigger.subplots_adjust(hspace=.0)


    # saving all the figures
    save_dir = '/'.join([
        analysis_dir,
        'analysis_plots',
        stream,
        simulation
    ])
    os.makedirs(save_dir, exist_ok=True)
    save_flag = True
    if save_flag:
        fig_temp.savefig(save_dir+'/fig_temp.png')
        fig_temp_heat.savefig(save_dir+'/fig_heat.png')
        fig_chi2.savefig(save_dir+'/fig_chi2_trig.png')
        fig_hist_trig.savefig(save_dir+'/fig_hist_trig.png')
        fig_hist_trig_ev.savefig(save_dir+'/fig_hist_trig_ev.png')
        fig_cross.savefig(save_dir+'/fig_cross.png')
        fig_ion.savefig(save_dir+'/fig_ion.png')
        fig_virtual.savefig(save_dir+'/fig_virtual.png')
        fig_10kev.savefig(save_dir+'/fig_10kev.png')
        fig_trigger.savefig(save_dir+'/fig_trigger.png')