#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 09:17:47 2019

@author: misiak
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import json
from plot_addon import basic_corner

# matplotlib configuration
plt.close('all')
plt.rcParams['text.usetex']=True

 

#stream = 'tg27l000'
save_flag = True

with open('/home/misiak/Data/data_run57_neutron/stream_config.json', 'r') as cfg:
    config = json.load(cfg)

calibration_quality = list()
for stream in config['Calibration']:
    save_dir = '/home/misiak/Analysis/fine_root_analysis/fond_neutron/{}'.format(stream)
    raw_load = np.load(save_dir+'/energy_quality.npy')
    data_quality = np.delete(raw_load, 1, axis=1)
    calibration_quality.append(data_quality)

energy_quality = np.concatenate(calibration_quality)

#energy_all = np.load(SAVE_DIR+'/energy_calibrated.npy')
#energy_quality = np.load(SAVE_DIR+'/energy_quality.npy')
#
## removing the heatB channel
#energy_all = np.delete(energy_all, 1, axis=1)
#energy_quality = np.delete(energy_quality, 1, axis=1)


heat, A, B, C, D = energy_quality.T

# =============================================================================
# FIDUCIAL CUT
# =============================================================================

fid_thresh = 0.5
fid_cut = (abs(A)<fid_thresh) & (abs(C)<fid_thresh)


fig, ax = plt.subplots()

ax.plot(A, C,
        label='quality events',
        ls='none', marker='.', color='k', alpha=0.3)

ax.plot(C[fid_cut], A[fid_cut],
        label='fiducial events',
        ls='none', marker='.', color='r', alpha=0.3)

ax.axhspan(-fid_thresh, fid_thresh,
           alpha=0.3, color='crimson', zorder=-5,
           label='fiducial cut')
ax.axvspan(-fid_thresh, fid_thresh,
           alpha=0.3, color='crimson', zorder=-5)

ax.set_xlim(-2, 12)
ax.set_ylim(-2, 12)

ax.set_ylabel('Ion C')
ax.set_xlabel('Ion A')
ax.grid()
ax.legend()


labels = ('heatA', 'ionA', 'ionB', 'ionC', 'ionD')
fig, axes = basic_corner(energy_quality, labels, color='k')
fig, axes = basic_corner(energy_quality[fid_cut], labels, axes=axes, color='r')


