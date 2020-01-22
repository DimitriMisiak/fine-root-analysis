#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 12:15:05 2020

@author: misiak
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from plot_addon import lighten_color

analysis_dir = '/home/misiak/Analysis/neutron_background'
analysis_simu_path = '/'.join([analysis_dir, 'simu_analysis.h5'])

df_simu = pd.read_hdf(
    analysis_simu_path,
    key='df',
    where='glitch_time_cut = True'
)

#%%
plt.close('all')

bins = np.arange(0, 50.5, 0.5)
bins_width = bins[1] - bins[0]
bins_array = bins[:-1] + (bins_width) / 2


# # =============================================================================
# # EFFICIENCY calculation
# # =============================================================================
energy_column = 'input_energy'

source_list = ['Calibration', 'Background']
simulation_list = ['flat_ER', 'flat_NR']

num_dict = dict()
eff_dict = dict()
for source in source_list:
    
    num_dict[source] = dict()
    eff_dict[source] = dict()
    for simulation in simulation_list:
        
        all_cut = (
            (df_simu['source'] == source)
            & (df_simu['simulation'] == simulation)
        )
        trigger_cut = df_simu['trigger_cut'] & all_cut
        quality_cut = df_simu['quality_cut'] & trigger_cut
        bulk_cut = df_simu['bulk_cut'] & quality_cut
       
        cut_dict = {
            'all': all_cut,
            'trigger': trigger_cut,
            'quality': quality_cut,
            'bulk': bulk_cut,
        }
        
        num_dict[source][simulation] = dict()
        eff_dict[source][simulation] = dict()
        for key, cut in cut_dict.items():
            num = np.histogram(
                df_simu[cut][energy_column],
                bins=bins
            )[0]
            num_dict[source][simulation][key] = num

            eff_dict[source][simulation][key] = (
                num / num_dict[source][simulation]['all']
            )


        
# ### PLOT number of events passing the cuts
# fig, axes = plt.subplots(nrows=2, ncols=2, num='Histogramm events passsing cuts',
#                          figsize=(10, 7), sharex='col', sharey='row')

# for im, source in enumerate(source_list):

#     ax = axes[im]

#     for js, simulation in enumerate(simulation_list):

#             a = ax[js]
            
#             for key, num in num_dict[source][simulation].items():
            
#                 line, = a.plot(
#                     bins_array,
#                     num,
#                     ls='steps-mid',
#                     marker='.',
#                     label=key
#                 )
#                 a.fill_between(
#                     bins_array,
#                     num,
#                     color=lighten_color(line.get_color()),
#                     step='mid',
#                 )
                
#             msg = '{} {} Events'.format(source, simulation).replace('_', ' ')
#             a.text(
#                 0.5, 0.1,
#                 msg,
#                 horizontalalignment='center',
#                 verticalalignment='center',
#                 transform=a.transAxes
#             )
#             a.grid()
#             a.set_yscale('log')

# a.legend(
#     loc='center left',
#     bbox_to_anchor=(1.05, 1.05)
# )
# axes[0,0].set_ylabel('Counts')
# axes[1,0].set_ylabel('Counts')
# axes[1,0].set_xlabel('Recoil Energy [keV]')
# axes[1,1].set_xlabel('Recoil Energy [keV]')

# fig.tight_layout()
# fig.subplots_adjust(hspace=0, wspace=0)


# ### PLOT number of events passing the cuts
# fig, axes = plt.subplots(nrows=2, ncols=2, num='Histogramm cut efficiency',
#                          figsize=(10, 7), sharex='col', sharey='row')

# for im, source in enumerate(source_list):

#     ax = axes[im]

#     for js, simulation in enumerate(simulation_list):

#             a = ax[js]
            
#             for key, eff in eff_dict[source][simulation].items():
            
#                 line, = a.plot(
#                     bins_array,
#                     eff,
#                     ls='steps-mid',
#                     marker='.',
#                     label=key
#                 )
#                 a.fill_between(
#                     bins_array,
#                     eff,
#                     color=lighten_color(line.get_color()),
#                     step='mid',
#                 )
                
#             msg = '{} {} Efficiency'.format(source, simulation).replace('_', ' ')
#             a.text(
#                 0.5, 0.1,
#                 msg,
#                 horizontalalignment='center',
#                 verticalalignment='center',
#                 transform=a.transAxes
#             )
#             a.grid()
#             a.set_yscale('log')

# a.legend(
#     loc='center left',
#     bbox_to_anchor=(1.05, 1.05)
# )
# axes[0,0].set_ylabel('Efficiency')
# axes[1,0].set_ylabel('Efficiency')
# axes[1,0].set_xlabel('Recoil Energy [keV]')
# axes[1,1].set_xlabel('Recoil Energy [keV]')

# fig.tight_layout()
# fig.subplots_adjust(hspace=0, wspace=0)

#%%
# =============================================================================
# CONTAMINATION calculation
# =============================================================================
energy_column = 'recoil_energy_bulk'

source_list = ['Calibration', 'Background']
simulation_list = ['flat_ER', 'line_1keV', 'line_10keV', 'flat_NR']

fine_cut = df_simu['trigger_cut'] & df_simu['quality_cut'] & df_simu['bulk_cut']
df_fine = df_simu[fine_cut]


gamma_cut = df_fine['gamma_cut']
neutron_cut = df_fine['neutron_cut']
ho_cut = df_fine['HO_cut']
all_cut = pd.Series(True, index=df_fine.index)

other_cut = ~(gamma_cut | neutron_cut | ho_cut)
pure_gamma_cut = gamma_cut & ~(neutron_cut | ho_cut)
pure_neutron_cut = neutron_cut & ~(gamma_cut | ho_cut)
pure_ho_cut = ho_cut & ~(gamma_cut | neutron_cut)

mix_gamma_neutron_cut = gamma_cut & neutron_cut & ~ho_cut
mix_gamma_ho_cut = gamma_cut & ho_cut & ~neutron_cut
mix_neutron_ho_cut = neutron_cut & ho_cut & ~gamma_cut

mix_gamma_neutron_ho_cut = gamma_cut & neutron_cut & ho_cut

pop_cut_dict = {
    'all': all_cut,
    'others': other_cut,
    'pure gamma': pure_gamma_cut,
    'pure neutron': pure_neutron_cut,
    'pure heat-only': pure_ho_cut,
    'mix gamma-neutron': mix_gamma_neutron_cut,
    'mix gamma-ho': mix_gamma_ho_cut,
    'mix neutron-ho': mix_neutron_ho_cut,
    'mix gamma-neutron-ho': mix_gamma_neutron_ho_cut
}


efficiency_link_dict = {
    'flat_ER': 'flat_ER',
    'line_1keV': 'flat_ER',
    'line_10keV': 'flat_ER',
    'flat_NR': 'flat_NR',
}

pop_dict = dict()
pop_corr_dict = dict()
for source in source_list:
    
    pop_dict[source] = dict()
    pop_corr_dict[source] = dict()
    for simulation in simulation_list:
        
        all_cut = (
            (df_simu['source'] == source)
            & (df_simu['simulation'] == simulation)
        )

        pop_dict[source][simulation] = dict()
        pop_corr_dict[source][simulation] = dict()
        for key, cut in pop_cut_dict.items():
            num = np.histogram(
                df_simu[all_cut & cut][energy_column],
                bins=bins
            )[0]
            
            pop_dict[source][simulation][key] = num
            key_eff = efficiency_link_dict[simulation]
            efficiency = eff_dict[source][key_eff]['bulk']
            pop_corr_dict[source][simulation][key] = num / efficiency
            

#%%
### PLOT number of events passing the cuts
fig, axes = plt.subplots(nrows=4, ncols=2, num='Stacked hist events population',
                          figsize=(10, 7), sharex='col', sharey='row')

###
for im, simulation in enumerate(simulation_list):

    ax = axes[im]
    for js, source in enumerate(source_list):

            a = ax[js]

            bot = bins_array * 0
            for key, num in pop_corr_dict[source][simulation].items():
                
                if key=='all':
                    a.plot(
                        bins_array,
                        num,
                        color='k',
                        ls='steps-mid',
                        zorder=10,
                        label=key
                    )
                    continue
                a.bar(
                    bins_array,
                    num,
                    width=bins_width,
                    bottom=bot,
                    label=key
                )
                bot += num

            msg = '{} {} Efficiency'.format(source, simulation).replace('_', ' ')
            a.text(
                0.5, 0.1,
                msg,
                horizontalalignment='center',
                verticalalignment='center',
                transform=a.transAxes
            )
            a.grid()

a.legend(
    loc='center left',
    bbox_to_anchor=(1.05, 1.05)
)

for ax in axes[:, 0]:
    ax.set_ylabel('Counts')
for ax in axes[-1, :]:
    ax.set_xlabel('Recoil Energy [keV]')

fig.tight_layout()
fig.subplots_adjust(hspace=0, wspace=0)







