#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 14:43:36 2020

@author: misiak
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.close('all')
plt.rcParams['text.usetex']=True
from tqdm import tqdm


analysis_dir = '/home/misiak/Analysis/neutron_background'
output_dir = '/'.join([analysis_dir, 'analysis_plots'])
extension='png'

h5type_list = [
    'data',
    'noise',
    'simu'
]

stream_list = [
    'tg18l005',
    'tg27l000',
    'tg28l000',
    'tg17l007',
    'tg19l010',
    'tg20l000',
    'tg21l000'
]

stream = stream_list[0]
h5type = 'data'
h5_path = '/'.join([ analysis_dir, "data_calibrated.h5" ])

df_analysis = pd.read_hdf(
    h5_path,
    key='df',
    where=(
        'stream = "{0}"'
    ).format(stream)
)

source = df_analysis['source'].unique()[0]
title = (
    '{0} {1} {2}'
).format(stream, h5type, source).replace('_', ' ') 


quality_cut = df_analysis['quality_cut']
A = df_analysis[quality_cut]['energy_adu_heat']
B = (1000 < A) & (A < 1400)

inf, med, sup = np.quantile(A[B], [0.16, 0.5, 0.84])

C = (inf < A) & (A < sup)

# D = df_analysis[quality_cut & C]
D = df_analysis[quality_cut]

D_a = D['energy_adu_ionA']
D_b = D['energy_adu_ionB']
D_c = D['energy_adu_ionC']
D_d = D['energy_adu_ionD']

top = (D_a - D_b) / (D_c + D_d)
bot = (D_c - D_d) / (D_a + D_b)

plt.figure()

plt.plot(
    bot,
    top,
    ls='none',
    marker='.'
)

plt.grid()
plt.axvline(-1, color='k')
plt.axvline(+1, color='k')
plt.axhline(-1, color='k')
plt.axhline(+1, color='k')