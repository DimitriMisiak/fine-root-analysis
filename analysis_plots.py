#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:04:43 2020

@author: misiak
"""

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
analysis_data_path = '/'.join([analysis_dir, 'data_analysis.h5'])

# stream_list = pd.read_hdf(
#     fine_data_path,
#     key='df',
#     columns=['stream',]
# )['stream'].unique()    

# from tqdm import tqdm
# for stream in stream_list:
if True:
    stream = 'tg28l000'

    df_analysis = pd.read_hdf(
        analysis_data_path,
        key='df',
        where='stream = "{}"'.format(stream)
    )
    
    source = df_analysis['source'].unique()[0]    
    
    
    # fig_temp = temporal_plot(stream, df_analysis)
    # fig_temp_heat = temporal_plot_heat_only(stream, df_analysis)
    
    
    # fig_chi2 = plot_chi2_vs_energy(stream, df_analysis)
    # # plotting the quality cuts
    # x_data = 10**np.linspace(-2, 5, int(1e4))
    # cut_ion = ion_chi2_threshold_function(
    #     analysis_parameters['ion_chi2_threshold'],
    #     x_data
    # )
    # cut_heat = heat_chi2_threshold_function(
    #     analysis_parameters[stream]['heat_chi2_threshold'],
    #     x_data)    
    # for i, ax in enumerate(fig_chi2.get_axes()):
    #     if i == 2:
    #         ax.plot(x_data, cut_heat, lw=1, color='k', label='quality cut')
    #     else:
    #         ax.plot(x_data, cut_ion, lw=1, color='k', label='quality cut')
    #     ax.set_xlim(10**-2, 10**5)
    #     ax.set_ylim(10**1, 10**9)
    #     ax.legend()            


    # fig_hist_trig = histogram_adu(stream, df_analysis)
    # # resize the plots
    # fig_hist_trig.get_axes()[0].set_xlim(-200, 2000)
    # for i, ax in enumerate(fig_hist_trig.get_axes()[:5]):
    #     if i==0:
    #         ax.set_xlim(-200, 2000)
    #     else:
    #         ax.set_xlim(-70, 70)    

  
    # samples = df_analysis[df_analysis.quality_cut][[
    #     'energy_adu_ionA',
    #     'energy_adu_ionB',
    #     'energy_adu_ionC',
    #     'energy_adu_ionD',
    # ]]
    # samples_corr = df_analysis[df_analysis.quality_cut][[
    #     'energy_adu_corr_ionA',
    #     'energy_adu_corr_ionB',
    #     'energy_adu_corr_ionC',
    #     'energy_adu_corr_ionD',
    # ]]
    # fig_cross, axes = basic_corner(
    #     samples.values,
    #     samples.columns,
    #     num = '{} ({}): Cross-talk Correction'.format(stream, source),
    #     label='raw',
    #     alpha=0.1,
    # )
    # basic_corner(
    #     samples_corr.values,
    #     samples_corr.columns,
    #     axes=axes,
    #     color='slateblue',
    #     zorder=-1,
    #     label='corrected'
    # )
    # for ax in fig_cross.get_axes():
    #     ax.axvline(0, color='r', zorder=-5)
    #     ax.axhline(0, color='r', zorder=-5)

    
    
    # fig_hist_trig_ev = histogram_ev(stream, df_analysis)
    # fig_ion = ion_vs_ion(stream, df_analysis)
    # fig_virtual = virtual_vs_virtual_ev(stream, df_analysis)
    
    # if save_flag:
    #     fig_temp.savefig(save_dir+'/fig_temp.png')
    #     fig_heat.savefig(save_dir+'/fig_heat.png')
        # fig_chi2.savefig(save_dir+'/fig_chi2_trig.png')
        # fig_hist_trig.savefig(save_dir+'/fig_hist_trig.png')
        # fig_cross.savefig(save_dir+'/fig_cross.png')
    