#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 12:18:01 2020

Mimic "ultimate_analysis.py".

@author: misiak
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from representation_functions import (
    band_cut_plots
)

analysis_dir = '/home/misiak/Analysis/neutron_background'

analysis_data_path = '/'.join([analysis_dir, 'data_analysis.h5'])
analysis_simu_path = '/'.join([analysis_dir, 'simu_analysis.h5'])


# df_data = pd.read_hdf(
#     analysis_data_path,
#     key='df',
# )

df_simu = pd.read_hdf(
    analysis_simu_path,
    key='df',
    where='glitch_time_cut = True'
)

#%%
plt.close('all')

# ### line 10keV
# bins = np.arange(5, 20, 0.1)

# ### line 1ekV
# bins = np.arange(0, 4, 0.01)

### flat
bins = np.arange(0, 51, 0.5)

bins_width = bins[1] - bins[0]
bins_array = bins[:-1] + (bins_width) / 2

all_cut = pd.Series(True, index=df_simu.index)
trigger_cut = df_simu['trigger_cut'] & all_cut
quality_cut = df_simu['quality_cut'] & trigger_cut
bulk_cut = df_simu['bulk_cut'] & quality_cut


plot_cut = (
    (df_simu.source == 'Calibration')
    & (df_simu.simulation == 'flat_NR')
)

df_pop = df_simu[bulk_cut & plot_cut]

gamma_cut = df_pop['gamma_cut']
neutron_cut = df_pop['neutron_cut']
ho_cut = df_pop['HO_cut']

all_cut = pd.Series(True, index=df_pop.index)

other_cut = ~(gamma_cut | neutron_cut | ho_cut)
pure_gamma_cut = gamma_cut & ~(neutron_cut | ho_cut)
pure_neutron_cut = neutron_cut & ~(gamma_cut | ho_cut)
pure_ho_cut = ho_cut & ~(gamma_cut | neutron_cut)

mix_gamma_neutron_cut = gamma_cut & neutron_cut & ~ho_cut
mix_gamma_ho_cut = gamma_cut & ho_cut & ~neutron_cut
mix_neutron_ho_cut = neutron_cut & ho_cut & ~gamma_cut

mix_gamma_neutron_ho_cut = gamma_cut & neutron_cut & ho_cut

pop_dict = {
    'others': other_cut,
    'pure gamma': pure_gamma_cut,
    'pure neutron': pure_neutron_cut,
    'pure heat-only': pure_ho_cut,
    'mix gamma-neutron': mix_gamma_neutron_cut,
    'mix gamma-ho': mix_gamma_ho_cut,
    'mix neutron-ho': mix_neutron_ho_cut,
    'mix gamma-neutron-ho': mix_gamma_neutron_ho_cut
}

quantity = 'recoil_energy_bulk'

x = list()
label = list()
for lab, cut in pop_dict.items():
    x.append( df_pop[cut][quantity] )
    label.append(lab)

fig, ax = plt.subplots()
ax.hist(
        x,
        bins=bins_array,
        label=label,
        histtype='stepfilled',
        stacked=True,
)

ax.hist(
        df_pop[all_cut][quantity],
        bins=bins_array,
        label='All',
        color='k',
        fill=False,
        histtype='step',
)

ax.legend()
ax.grid()


band_cut_plots('hey', df_pop)

# fig, ax = plt.subplots()

# ax.hist(
#         df_simu[plot_cut & bulk_cut]['input_energy'],
#         bins=bins_array,
#         alpha=0.3,
#         label='input energy'
# )

# ax.hist(
#         df_simu[plot_cut & bulk_cut]['recoil_energy_bulk'],
#         bins=bins_array,
#         alpha=0.3,
#         label='recoil energy',
# )

# ax.legend()
# ax.grid()


# def plot_eff(ax, bins, source, simulation, band_cut,
#              quantity='input_energy'):

#         plot_cut = (
#             (df_simu.source == source)
#             & (df_simu.simulation == simulation)
#         )
    
#         n_all, _, __ = ax.hist(
#             df_simu[plot_cut & all_cut][quantity],
#             bins=bins,
#             label='All events'
#         )

#         n_trigger, _, __ = ax.hist(
#             df_simu[plot_cut & trigger_cut][quantity],
#             bins=bins,
#             label='+ trigger'
#         )

#         n_quality, _, __ = ax.hist(
#             df_simu[plot_cut & quality_cut][quantity],
#             bins=bins,
#             label='+ quality'
#         )

#         n_bulk, _, __ = ax.hist(
#             df_simu[plot_cut & bulk_cut][quantity],
#             bins=bins,
#             label='+ bulk'
#         )

#         band_cut_truth = df_simu[band_cut]

#         n_band, _, __ = ax.hist(
#             df_simu[plot_cut & band_cut_truth & bulk_cut][quantity],
#             bins=bins,
#             label='+ {}'.format(band_cut).replace('_', ' ')
#         )


#         ax.legend(
#             loc='upper right',
#             title='{} SIMU {}'.format(source, simulation).replace('_', ' ')
#         )
#         ax.grid()
#         ax.set_yscale('log')
        
#         return None

# fig, ax = plt.subplots()

# plot_eff(ax, bins_array, 'Background', 'flat_ER', 'gamma_cut',
#           quantity='recoil_energy_bulk')


# =============================================================================
# QUENCHING PLOT WITH BACKGROUND AND CALIBRATION COMPARISON
# =============================================================================
# fine_data_cut = (
#     df_data['glitch_time_cut']
#     & df_data['energy_cut']
#     & df_data['quality_cut']
#     & df_data['bulk_cut']
# )

# df_calibration = df_data[fine_data_cut & (df_data.source == 'Calibration')]
# df_background = df_data[fine_data_cut & (df_data.source == 'Background')]

# fig_quenching, ax = plt.subplots(
#     figsize=(10, 7),
#     num='Quenching plot comparison'
# )

# ax.plot(
#         df_calibration['recoil_energy_bulk'],
#         df_calibration['quenching_bulk'],
#         ls='none', marker='.',
#         color='deepskyblue',
#         label='Calibration',
#         alpha=0.5,
# )

# ax.plot(
#         df_background['recoil_energy_bulk'],
#         df_background['quenching_bulk'],
#         ls='none', marker='.',
#         color='k',
#         label='Background',
#         alpha=0.5,
# )

# ax.legend()
# ax.set_xlabel('Recoil Energy [keV]')
# ax.set_ylabel('Quenching Factor [fraction]')
# ax.set_ylim(-0.25, 1.5)
# ax.set_xlim(0, 30)
# ax.grid(True)
# fig_quenching.tight_layout()

# =============================================================================
# CUT EFFICIENCY
# =============================================================================

# all_cut = df_simu['glitch_time_cut']# & df_simu['energy_cut']
# all_cut = pd.Series(True, index=df_simu.index)
# trigger_cut = df_simu['trigger_cut'] & all_cut
# quality_cut = df_simu['quality_cut'] & trigger_cut
# bulk_cut = df_simu['bulk_cut'] & quality_cut

# source_configurations = ['Calibration', 'Background']

# # simulation_list = ['flat_ER', 'flat_NR', 'line_1keV', 'line_10keV']
# simulation_list = ['flat_ER', 'flat_NR']

# def plot_eff(ax, bins, source, simulation):

#         plot_cut = (
#             (df_simu.source == source)
#             & (df_simu.simulation == simulation)
#         )
    
#         n_all, _, __ = ax.hist(
#             df_simu[plot_cut & all_cut]['input_energy'],
#             bins=bins,
#             label='All events'
#         )

#         n_trigger, _, __ = ax.hist(
#             df_simu[plot_cut & trigger_cut]['input_energy'],
#             bins=bins,
#             label='+ trigger'
#         )

#         n_quality, _, __ = ax.hist(
#             df_simu[plot_cut & quality_cut]['input_energy'],
#             bins=bins,
#             label='+ quality'
#         )

#         n_bulk, _, __ = ax.hist(
#             df_simu[plot_cut & bulk_cut]['input_energy'],
#             bins=bins,
#             label='+ bulk'
#         )

#         if simulation == 'flat_ER':
#             band_cut = df_simu['gamma_cut']
#         elif simulation == 'flat_NR':
#             band_cut = df_simu['neutron_cut']

#         n_band, _, __ = ax.hist(
#             df_simu[plot_cut & band_cut & bulk_cut]['input_energy'],
#             bins=bins,
#             label='+ {}'.format(simulation).replace('_', ' ')
#         )


#         a.legend(loc='upper right', title='{} SIMU {}'.format(mode, stype).replace('_', ' '))
#         a.grid()
#         a.set_yscale('log')

#         eff_local = dict()
#         eff_local['trigger'] = n_trigger / n_all
#         eff_local['quality'] = n_quality / n_all
#         eff_local['bulk'] = n_bulk / n_all
#         eff_local[simulation] = n_band / n_all
#         # eff_dict[source][stype] = eff_local
        
#         return None


# bins = np.arange(0, 51, 1)
# bins_width = bins[1] - bins[0]
# bins_array = bins[:-1] + (bins_width) / 2
# eff_x_array = bins_array

# eff_dict = dict()

# ###
# fig, axes = plt.subplots(nrows=2, ncols=2, num='Hist E input',
#                          figsize=(10, 7), sharex='col', sharey='row')

# for im, mode in enumerate(['Background', 'Calibration']):

#     ax = axes[im]
#     eff_dict[mode] = dict()

#     for js, stype in enumerate(simulation_list):

#         a = ax[js]
#         plot_eff(a, bins_array, mode, stype)


# axes[0,0].set_ylabel('Counts')
# axes[1,0].set_ylabel('Counts')
# axes[1,0].set_xlabel('Input Energy [keV]')
# axes[1,1].set_xlabel('Input Energy [keV]')

# fig.tight_layout()
# fig.subplots_adjust(hspace=0, wspace=0)