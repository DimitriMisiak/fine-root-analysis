#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 13:51:17 2019

@author: misiak
"""

import numpy as np
import matplotlib.pyplot as plt

all_data = np.load(
        '/home/misiak/Analysis/fine_root_analysis/fond_neutron/tg17l007/Flat_Analytical_SimuOnly_0.0000_50.0000_ER/input_energy_all.npy'
)

trigger_data = np.load('/home/misiak/Analysis/fine_root_analysis/fond_neutron/tg17l007/Flat_Analytical_SimuOnly_0.0000_50.0000_ER/input_energy_trigger.npy')

quality_data = np.load('/home/misiak/Analysis/fine_root_analysis/fond_neutron/tg17l007/Flat_Analytical_SimuOnly_0.0000_50.0000_ER/input_energy_quality.npy')


plt.close('all')

n, bins, patches = plt.hist(all_data[:], bins=1000, label='All data', alpha=0.3)

plt.hist(trigger_data[:], bins=bins, label='Trigger Data', alpha=0.3)

plt.hist(quality_data[:], bins=bins, label='Quality Data', alpha=0.3)

plt.legend()