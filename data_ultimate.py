#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 13:34:35 2020

@author: misiak
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

plt.close('all')

analysis_dir = '/home/misiak/Analysis/neutron_background'
simu_path = '/'.join([analysis_dir, 'simu_analysis.h5'])
data_path = '/'.join([analysis_dir, 'data_analysis.h5'])


df_data = pd.read_hdf(
    data_path,
    key='df',
    where='source = "Calibration" & quality_cut = True'
)

df_simu = pd.read_hdf(
    simu_path,
    key='df',
    where=(
        'source = "Calibration"'
        '& simulation = "line_1keV"'
        '& quality_cut = True'
        '& bulk_cut = True'
        '& energy_cut = True'
        '& trigger_cut = True'
    )
)

neutron_cut = df_simu['neutron_cut']
gamma_cut = df_simu['gamma_cut']
HO_cut = df_simu['HO_cut']
bonus_cut = neutron_cut & gamma_cut

fig, ax = plt.subplots()
for suffix in ['neutron_cut', 'gamma_cut', 'bonus_cut']:
    
    if suffix == 'neutron_cut':
        color='blue'
        cut = df_simu[suffix]
    if suffix == 'gamma_cut':
        color='orange'
        cut = df_simu[suffix]
    if suffix == 'bonus_cut':
        color='k'
        cut = bonus_cut
        
    label = suffix.replace('_', ' ')
    ax.plot(
            df_simu[cut]['energy_heat'],
            df_simu[cut]['energy_ion_total'],
            ls='none',
            marker='.',
            label=label,
            alpha=0.5,
            color=color
    )
ax.grid()
ax.legend()


# quality_cut = df_data['quality_cut']
# bulk_cut = df_data['bulk_cut']
# energy_cut = df_data['energy_cut']
# fine_cut = quality_cut & bulk_cut & energy_cut

# neutron_cut = df_data['neutron_cut']
# gamma_cut = df_data['gamma_cut']
# HO_cut = df_data['HO_cut']

# fig, ax = plt.subplots()
# for suffix in ['neutron_cut', 'gamma_cut', 'HO_cut']:
    
#     cut = df_data[suffix]
#     label = suffix.replace('_', ' ')
#     ax.plot(
#             df_data[fine_cut][cut]['energy_heat'],
#             df_data[fine_cut][cut]['energy_ion_total'],
#             ls='none',
#             marker='.',
#             label=label
#     )
# ax.legend()

# fig, ax = plt.subplots()
# bins = np.arange(0, 51, 1)
# ax.hist(
#     df_data[fine_cut & neutron_cut]['recoil_energy_bulk'],
#     bins=bins
# )
# ax.set_yscale('log')


# # =============================================================================
# # EFFICIENCY
# # plt.close('all')
# # =============================================================================

# mode = 'Background'
# stype = 'NR'
# #bins = 50
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

#     for js, stype in enumerate(['ER', 'NR']):

#         a = ax[js]

#         n_all, _, __ = a.hist(D[mode]['simu']['input_all'][stype], bins=bins, label='all')

#         n_trig,_, __ = a.hist(D[mode]['simu']['input_trigger'][stype], bins=bins, label='trigger')

#         n_qual,_, __ = a.hist(D[mode]['simu']['input_quality'][stype], bins=bins, label='quality')

#         n_fid,_, __ = a.hist(D[mode]['simu']['input_fid'][stype], bins=bins, label='fid')

#         if stype == 'NR':
#             band = 'input_neutron'
#         elif stype == 'ER':
#             band = 'input_gamma'
#         n_band, _, __ = a.hist(D[mode]['simu'][band][stype], bins=bins, label='NR')

#         a.legend(loc='upper right', title='{} SIMU {}'.format(mode, stype))
#         a.grid()
#         a.set_yscale('log')

#         eff_local = dict()
#         eff_local['trigger'] = n_trig / n_all
#         eff_local['quality'] = n_qual / n_all
#         eff_local['fid'] = n_fid / n_all
#         eff_local['band'] = n_band / n_all
#         eff_dict[mode][stype] = eff_local

# axes[0,0].set_ylabel('Counts')
# axes[1,0].set_ylabel('Counts')
# axes[1,0].set_xlabel('Input Energy [keV]')
# axes[1,1].set_xlabel('Input Energy [keV]')

# fig.tight_layout()
# fig.subplots_adjust(hspace=0, wspace=0)

# ###
# fig, axes = plt.subplots(nrows=2, ncols=2, num='Efficiency as f(E input)',
#                          figsize=(10, 7), sharex='col', sharey='row')

# for im, mode in enumerate(['Background', 'Calibration']):

#     ax = axes[im]

#     for js, stype in enumerate(['ER', 'NR']):

#         a = ax[js]

#         for eff_key, eff_array in eff_dict[mode][stype].items():
#             line, = a.plot(
#                     bins_array, eff_array,
#                     ls='steps-mid', alpha=0.3,
#             )
#             a.errorbar(
#                     bins_array, eff_array,
#                     xerr = bins_width/2, yerr = eff_array*0.05,
#                     label=eff_key,
#                     ls='none', marker='.',
#                     color=line.get_color()
#             )

#         a.legend(loc='upper right', title='{} SIMU {}'.format(mode, stype))
#         a.grid()
#         a.set_yscale('log')

# axes[0,0].set_ylabel('Efficiency Fraction')
# axes[1,0].set_ylabel('Efficiency Fraction')
# axes[1,0].set_xlabel('Input Energy [keV]')
# axes[1,1].set_xlabel('Input Energy [keV]')

# fig.tight_layout()
# fig.subplots_adjust(hspace=0, wspace=0)


# #%%
# # =============================================================================
# # BACKGROUND
# #plt.close('all')
# # =============================================================================

# mode = 'Background'
# stype = 'NR'
# bins = np.arange(0, 51, 0.5)
# bins_width = bins[1] - bins[0]
# bins_array = bins[:-1] + (bins_width) / 2

# adv_eff_dict = dict()
# for mode, eff_mode in eff_dict.items():
#     adv_eff_dict[mode] = dict()

#     for stype, eff_stype in eff_mode.items():
#         adv_eff_dict[mode][stype] = dict()

#         for ctype, eff_ctype in eff_stype.items():
#             adv_eff_dict[mode][stype][ctype] = np.interp(
#                     bins_array,
#                     eff_x_array,
#                     eff_ctype
#             )


# ###
# fig, axes = plt.subplots(nrows=2, ncols=2, num='Energy Spectrum with Cuts',
#                          figsize=(10, 7), sharex='col', sharey='row')

# adv_data_dict = dict()

# for im, mode in enumerate(['Background', 'Calibration']):

#     ax = axes[im]
#     adv_data_dict[mode] = dict()

#     for js, stype in enumerate(['ER', 'NR']):

#         a = ax[js]

#         for data_key, data_array in D[mode]['data'].items():

#             if (data_key == 'neutron') and (stype == 'ER'):
#                 continue

#             if (data_key == 'gamma') and (stype == 'NR'):
#                 continue

#             n ,bins, patches = a.hist(
#                     gen_er(data_array),
#                     bins=bins,
#                     label=data_key
#             )

#             if (data_key == 'gamma') or (data_key == 'neutron'):
#                 adv_data_dict[mode][stype] = n

#         a.legend(loc='upper right', title='{} SIMU {}'.format(mode, stype))
#         a.grid()
#         a.set_yscale('log')

# axes[0,0].set_ylabel('Counts')
# axes[1,0].set_ylabel('Counts')
# axes[1,0].set_xlabel('Recoil Energy [keV]')
# axes[1,1].set_xlabel('Recoil Energy [keV]')

# fig.tight_layout()
# fig.subplots_adjust(hspace=0, wspace=0)


# ### NORMALIZATION
# ### COUNTS TO DRU

# mass_ge = 0.038 #kg

# exposure_dict = {
#         'Background':37.45/24, #days
#         'Calibration':76./24, #days
# }

# DRU_dict = dict()
# inte_dict = dict()

# fig, axes = plt.subplots(nrows=2, ncols=2, num='Energy Spectrum in DRU',
#                          figsize=(10, 8), sharex='col', sharey='row')

# for im, mode in enumerate(['Background', 'Calibration']):

#     ax = axes[im]
#     DRU_dict[mode] = dict()

#     for js, stype in enumerate(['ER', 'NR']):

#         a = ax[js]

#         data_bin_array = adv_data_dict[mode][stype]
#         eff_array = adv_eff_dict[mode][stype]['band']
#         exposure = exposure_dict[mode]

#         DRU_dict[mode][stype] = data_bin_array / (eff_array * exposure * bins_width * mass_ge)

#         a.plot(
#                 bins_array,
#                 DRU_dict[mode][stype],
#                 label='{} {} [DRU]'.format(mode, stype),
#                 ls='steps-mid',
#         )

#         if stype == 'NR':
#             cut_2kev = (bins_array >= 2)

#             inte = np.trapz(
#                     DRU_dict[mode][stype][cut_2kev],
#                     bins_array[cut_2kev]
#             )

#             inte_dict[mode] = inte

#             a.fill_between(
#                     bins_array[cut_2kev],
#                     DRU_dict[mode][stype][cut_2kev],
#                     label='In [2keV-50keV]:\n{:.1f} Counts/kg/days'.format(inte),
#                     step='mid',
#                     color='coral',

#             )

#         a.legend(loc='upper right', title='{} SIMU {}'.format(mode, stype))
#         a.grid()
#         a.set_yscale('log')

# axes[0,0].set_ylabel('DRU [Efficiency corrected Counts/keV/kg/days]')
# axes[1,0].set_ylabel('DRU [Efficiency corrected Counts/keV/kg/days]')
# axes[1,0].set_xlabel('Recoil Energy [keV]')
# axes[1,1].set_xlabel('Recoil Energy [keV]')


# for ax in axes:
#     for a in ax:
#         a.set_xscale('log')

# inte_bkgd = inte_dict['Background']
# inte_calib = inte_dict['Calibration']
# ratio = inte_calib / inte_bkgd

# fig.suptitle(
#         (
#                 'In ROI [2keV-50keV]: {:.1f} Calib / {:.1f} Bkgd = {:.2f} Ratio'
#         ).format(inte_calib, inte_bkgd, ratio)
# )

# fig.tight_layout(rect=(0, 0, 1, 0.950))
# fig.subplots_adjust(hspace=0, wspace=0)


# #%%
# ### Money plot
# calib_er = DRU_dict['Calibration']['ER']
# calib_nr = DRU_dict['Calibration']['NR']
# bkgd_er = DRU_dict['Background']['ER']
# bkgd_nr = DRU_dict['Background']['NR']

# calib_tot = calib_er + calib_nr
# bkgd_tot = bkgd_er + bkgd_nr

# array_list = [
#     # calib_tot,
#     # bkgd_tot,
#     calib_er,
#     calib_nr,
#     bkgd_er,
#     bkgd_nr
# ]

# color_list = [
#     # 'grey',
#     # 'lightgrey',
#     'red',
#     'blue',
#     'coral',
#     'deepskyblue',
#     ]

# legend_list =[
#     'Calibration ER band',
#     'Calibration NR band',
#     'Background ER band',
#     'Background NR band',
# ]

# fig, ax = plt.subplots()

# for i,dru_array in enumerate(array_list):

#     c = color_list[i]
#     leg = legend_list[i]
    
#     zorder=1
    
#     if i in (0,2):
#         zorder=5
    
#     ax.plot(
#         bins_array,
#         dru_array,
#         ls='steps-mid',
#         alpha=1,
#         color=c,
#         lw=3,
#         zorder=zorder,
#         #path_effects=cartoon,
#         label=leg
#     )

#     ax.plot(
#         bins_array,
#         dru_array,
#         ls='none',
#         alpha=1,
#         marker='.',
#         color=c,
#         zorder=zorder,
#         #path_effects=cartoon,
#     )

#     if i in (1,3):
        
#         if i == 1:
#             msg = (
#                 'Calibration events in [2keV-50keV]:\n{:.2e} Counts/kg/days'
#             ).format(inte_dict['Calibration'])
#         if i == 3:
#             msg = (
#                 'Background events in [2keV-50keV]:\n{:.2e} Counts/kg/days'
#             ).format(inte_dict['Background'])
            
#         ax.fill_between(
#             bins_array,
#             dru_array,
#             step='mid',
#             color=lighten_color(c),
#             zorder=-1,
#             label=msg
#         )

# ax.axvspan(0, 2, color='k', alpha=0.3, zorder=5)
# ax.axvline(2, ls='--', color='k', zorder=5, 
#            label='Analysis Threshold: 2keV')



# ax.legend(handler_map={str: LegendTitle()})

# ax.set_xlim(0.25, 50)
# ax.set_ylim(1e2, 1e7)
# ax.set_yscale('log')
# ax.set_ylabel('DRU [Efficiency corrected Counts/keV/kg/days]')
# ax.set_xlabel('Recoil Energy [keV]')

# ax.grid(which='both', alpha=0.5)
# fig.tight_layout()

