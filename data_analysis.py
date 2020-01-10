#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:55:53 2020

@author: misiak
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

analysis_dir = '/home/misiak/Analysis/neutron_background'
h5_path =  '/home/misiak/Analysis/neutron_background/data.h5'
simu_path =  '/home/misiak/Analysis/neutron_background/simu.h5'
fine_path = '/home/misiak/Analysis/neutron_background/data_useful.h5'

# =============================================================================
# FOR DEBUG, to be deleted later
# =============================================================================
# for the production

# df_iterator = pd.read_hdf(
#     h5_path,
#     key='data',
#     chunksize=50000
# )

# for i, df in enumerate(df_iterator):
#     pass

# for the moment, one small df to build the code:

df = pd.read_hdf(
    h5_path,
    key='data',
    start=0,
    stop=50000
)

# =============================================================================
# INITIALIZATION of the FINE DataFrame
# =============================================================================
df_fine = pd.DataFrame()
for col in ['detector', 'source', 'stream']:
    df_fine[col] = df[col]

# =============================================================================
# TIMESTAMP
# row by row
# columns needed: ['MicroStp_filt_decor', 'NumPart_filt_decor', 'frequency']
# columns returned: ['timestamp']
# =============================================================================
df_fine['timestamp'] = (
    df.MicroStp_filt_decor / (df.frequency * 3600)
    + df.NumPart_filt_decor
)

# =============================================================================
# CUSTOM TIME CUT
# row by row
# (stream by stream, graphical)
# columns needed: ['stream']
# columns returned: ['custom_time_cut']

# to be implemented after processing
# =============================================================================
# custom_time_cut_json_path = '/'.join([analysis_dir, 'custom_time_cut.json'])
# with open(custom_time_cut_json_path, "r") as cfg:
#     custom_time_cut_dict = json.load(cfg)

# if stream == 'tg18l005':
#     time_cut = (7.4, 7.6)

# if stream == 'tg27l000':
#     thresh_heat = 700
#     time_cut = (7, 11.3)

# if stream == 'tg28l000':
#     thresh_heat = 700
#     time_cut = (7.4, 8.05)

# custom_time_cut_series = pd.Series(
#     np.ones(df_temp.shape[0], dtype=bool),
#     index=df_temp.index
# )

# =============================================================================
# MAINTENANCE CUT
# row by row
# (stream by stream, graphical)
# columns needed: ['maintenance_duration', 'maintenance_cycle']
# columns returned: ['maintenance_cut']
# =============================================================================
full_maintenance_duration = (
    (df.maintenance_duration + df.maintenance_cycle)/3600
)
remainder = df_fine.timestamp % full_maintenance_duration

df_fine['maintenance_cut'] = (
    remainder > (df.maintenance_duration/3600)
)

# =============================================================================
# CHI2 CUT
# row by row
# (stream by stream, graphical)
# columns needed: ['stream', 'Chi2'+] + chosen cut
# columns returned: ['chi2_heat_cut', 'chi2_ion_cut']
# =============================================================================

# =============================================================================
# EXTRACTING CHI2
# =============================================================================
chi2_key_dict = {
    'chi2_heat': 'chi2_OF[0]_filt_decor',
    'chi2_ionA': 'chi2_OF[2]_filt_decor',
    'chi2_ionB': 'chi2_OF[3]_filt_decor',
    'chi2_ionC': 'chi2_OF[4]_filt_decor',
    'chi2_ionD': 'chi2_OF[5]_filt_decor',
}

for col, key in chi2_key_dict.items():
    df_fine[col] = df[key]

# =============================================================================
# ION OFFSET CUT
# row by row
# (stream by stream, graphical)
# columns needed: ['Offset'] + chosen cut
# columns returned: ['ion_offset_cut']
# =============================================================================

# =============================================================================
# EXTRACTING OFFSET
# =============================================================================
offset_key_dict = {
    'offset_heat': 'Off[0]_raw',
    'offset_ionA': 'Off[2]_raw',
    'offset_ionB': 'Off[3]_raw',
    'offset_ionC': 'Off[4]_raw',
    'offset_ionD': 'Off[5]_raw',
}

for col, key in offset_key_dict.items():
    df_fine[col] = df[key]

# =============================================================================
# EXTRACTING SLOPE
# =============================================================================
slope_key_dict = {
     'slope_ionA': 'Slope_Ion[0]_raw',
     'slope_ionB': 'Slope_Ion[1]_raw',
     'slope_ionC': 'Slope_Ion[2]_raw',
     'slope_ionD': 'Slope_Ion[3]_raw',
}

for col, key in offset_key_dict.items():
    df_fine[col] = df[key]

# =============================================================================
# QUALITY CUT
# row by row
# columns needed: [custom_time_cut, maintenance_cut, chi2_heat_cut,
#                  chi2_ion_cut, ion_offset_cut]
# columns returned: ['quality_cut']
# =============================================================================

# =============================================================================
# CROSSTALK CORRECTION
# row by row
# (all data, graphical)
# columns needed: [Energy_OF Ion] + crosstalk matrix
# columns returned: ['Energy_adu_crosstalk', 'Energy_adu']
# =============================================================================

# =============================================================================
# CALIBRATION
# row by row
# (stream by stream, graphical)
# columns needed: [Energy_OF] + position of the calibration peaks
# columns returned: ['Energy']
# =============================================================================

# =============================================================================
# VIRTUAL IONIZATION CHANNEL
# row by row
# (all data, graphical)
# columns needed: ['Energy', 'Polar']
# columns returned: ['Energy_Total', 'Charge_Conservation', 'Collect', 'Veto']
# =============================================================================

# =============================================================================
# EXTRACTING ENERGY
# =============================================================================
energy_key_dict = {
    'energy_adu_heat': 'Energy_OF[0]_filt_decor',
    'energy_adu_ionA': 'Energy_OF[2]_filt_decor',
    'energy_adu_ionB': 'Energy_OF[3]_filt_decor',
    'energy_adu_ionC': 'Energy_OF[4]_filt_decor',
    'energy_adu_ionD': 'Energy_OF[5]_filt_decor',
}

for col, key in energy_key_dict.items():
    df_fine[col] = df[key]

# =============================================================================
# FINALIZATION
# =============================================================================
df_fine.to_hdf(fine_path, key='df', mode='w', format='table')
