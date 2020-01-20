#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 14:30:17 2020

@author: misiak
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from data_analysis import (
    analysis_parameters,
    ion_channel_labels,
)

from plot_addon import basic_corner

stream = 'tg28l000'
source = 'Background'

analysis_dir = '/home/misiak/Analysis/neutron_background'
# hdf5_path = '/'.join([analysis_dir, 'data_analysis.h5'])
hdf5_path = '/'.join([analysis_dir, 'simu_analysis.h5'])

df = pd.read_hdf(
    hdf5_path,
    key='df',
    # where='stream = "{}"'.format(stream)
    where='source = "{}" & simulation = "line_1keV"'.format(source)
)

quality_cut = df['quality_cut']
bulk_cut = df['bulk_cut']
guard_cut = df['guard_cut']
energy_cut = df['energy_cut']
neutron_cut = df['neutron_cut']
gamma_cut = df['gamma_cut']
HO_cut = df['HO_cut']

# simu only
trigger_cut = df['trigger_cut']

plt.close('all')
# =============================================================================
# BULK, GUARD, all
# =============================================================================
# cut_dict = {
#     "all": energy_cut & quality_cut,
#     "bulk": energy_cut & quality_cut & bulk_cut,
#     "guard": energy_cut & quality_cut & guard_cut,
    
# }

# color_dict = {
#     "all": 'k',
#     "bulk": 'blue',
#     "guard": 'red',
# }    

# flag = 0
# for label, cut in cut_dict.items():
 
#     samples = df[cut][[
#         'energy_ionA',
#         'energy_ionB',
#         'energy_ionC',
#         'energy_ionD',
#     ]]
    
#     if flag == 0:
#         axes=None
        
#     fig_cross, axes = basic_corner(
#         samples.values,
#         samples.columns,
#         num = '{}: bulk, guard, all'.format(stream),
#         label=label,
#         axes=axes,
#         color=color_dict[label],
#     )
    
#     flag += 1

# =============================================================================
# NEUTRON, GAMMA, HO, ALL
# =============================================================================

all_cut = energy_cut & quality_cut & bulk_cut & trigger_cut
cut_dict = {
    "all": all_cut,
    "neutron": all_cut & neutron_cut,
    "gamma": all_cut & gamma_cut,
    "HO": all_cut & HO_cut,
}

# simu only
cut_dict = {k: v & trigger_cut for k,v in cut_dict.items()}

color_dict = {
    "all": 'k',
    "neutron": 'blue',
    "gamma": 'red',
    "HO": 'green'
}    

mks = 10
fig, ax = plt.subplots()
for label, cut in cut_dict.items():
    
    ei_array = df[cut]['energy_ion_bulk']
    ec_array = df[cut]['energy_heat']
    
    ax.plot(
        ec_array,
        ei_array,
        ls='none', marker='o',
        markersize=mks,
        color=color_dict[label],
        label=label
    )
    mks -= 2
ax.legend()
ax.grid()

neutron_cut = cut_dict['neutron']
gamma_cut = cut_dict['gamma']
HO_cut = cut_dict['HO']

pure_gamma = gamma_cut & ~(neutron_cut)
pure_neutron = neutron_cut & ~(gamma_cut)
mix_gamma_neutron = neutron_cut & gamma_cut
not_gamma_neutron = all_cut & ~(neutron_cut | gamma_cut)

print('Population counts:')
print('------------------')

print('All events = {}'.format(all_cut.sum()))

print('Pure gamma = {}'.format( pure_gamma.sum() ))
print('Pure neutron = {}'.format( pure_neutron.sum() ))
print('Mix gamma-neutron = {}'.format( mix_gamma_neutron.sum() ))
print('Not gamma-neutron = {}'.format( not_gamma_neutron.sum() ))
print('Sum 4 previous = {}'.format(
    (pure_gamma | pure_neutron | mix_gamma_neutron | not_gamma_neutron).sum()
))
# print('Pure gamma = {}'.format( ( gamma_cut & ~(neutron_cut) ).sum() ))
# print('Pure neutron = {}'.format( ( neutron_cut & ~(gamma_cut) ).sum() ))
# print('Mix gamma-neutron = {}'.format( (neutron_cut & gamma_cut).sum() ))
# print('Not gamma-neutron = {}'.format( (all_cut & ~(neutron_cut | gamma_cut)).sum() ))


