#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:22:18 2020

@author: misiak
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from plot_addon import lighten_color, LegendTitle

analysis_dir = '/home/misiak/Analysis/neutron_background'
analysis_simu_path = '/'.join([analysis_dir, 'simu_science.h5'])
analysis_data_path = '/'.join([analysis_dir, 'data_science.h5'])

# source_list = ['Background', 'Calibration']

df_data = pd.read_hdf(
    analysis_data_path,
    key='df',
    where=(
        'glitch_time_cut == True'
        '& quality_cut == True'
        '& bulk_cut == True'
        '& charge_conservation_cut == True'
    )
)

#%%
plt.close('all')

stream_list = df_data['stream'].unique()
stream_list.sort()

ho_cut = df_data['HO_cut']
df_ho = df_data[ho_cut]

fig, axes = plt.subplots(
    nrows=len(stream_list),
)

for ax, stream in zip(axes, stream_list):
    df = df_ho[df_ho['stream'] == stream]
    
    ax.hist(
            df['timestamp'],
            bins=100
    )
    
    ax.set_ylabel(stream)