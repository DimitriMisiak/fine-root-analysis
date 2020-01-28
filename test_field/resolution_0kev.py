#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:46:53 2020

@author: misiak
"""


import pandas as pd


analysis_noise_path = '/home/misiak/Analysis/neutron_background/noise_analysis.h5'

df_noise = pd.read_hdf(analysis_noise_path)

channel_list = [
    'energy_heat',
    'energy_ionA',
    'energy_ionB',
    'energy_ionC',
    'energy_ionD',
    'energy_ion_bulk',
    'energy_ion_total'
]

source_list = ['Background', 'Calibration']

for source in source_list:
    
    df = df_noise[(df_noise['source'] == source) & df_noise['quality_cut']]
    
    print('for {}'.format(source))
    
    for channel in channel_list:
        
        std = df[channel].std()
        
        print('{}: {}'.format(channel, std))