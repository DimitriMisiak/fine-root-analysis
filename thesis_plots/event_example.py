#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 11:03:03 2020

@author: misiak
"""


import numpy as np
import matplotlib.pyplot as plt
import uproot
import pandas as pd

stream = "tg27l000"

# h5_path = "/home/misiak/Analysis/neutron_background/data_science.h5"

# df = pd.read_hdf(h5_path)

plt.close('all')
plt.rcParams['text.usetex']=True

data_path = "/home/misiak/Data/data_run57_neutron/Data/tg27l000/RED80/RootFiles/TriggerData_tg27l000_S04_RED80_ChanTrig0.root"

tree = uproot.open(data_path)['tree']

heat = tree['Trace_Heat_A_Raw'].array()
ionA = tree['Trace_Ion_A_Raw'].array()
ionB = tree['Trace_Ion_B_Raw'].array()
ionC = tree['Trace_Ion_C_Raw'].array()
ionD = tree['Trace_Ion_D_Raw'].array()

channels = [heat, ionA, ionB, ionC, ionD]

time_array = np.arange(0, 1, 200**-1)

ndim = 5
fig, axes = plt.subplots(nrows=ndim, figsize=(10, 8))

ind = 10
cmap = plt.get_cmap('jet')

# for i in range(ind):
# for i in [np.random.randint(1, 7000),]:
for i in [3265,]:

    for n in range(ndim):
        
        ax = axes[n]
        chan = channels[n]
        
        ax.plot(
            time_array,
            chan[i],
            color=cmap(i/ind)
        )
        
for ax in axes:
    # ax.set_ylim(-32000, 32000)
    ax.grid()

axes[0].set_ylabel('Heat', labelpad=10)
axes[1].set_ylabel('Ion. A', labelpad=10)
axes[2].set_ylabel('Ion. B', labelpad=10)
axes[3].set_ylabel('Ion. C', labelpad=10)
axes[4].set_ylabel('Ion. D', labelpad=10)
axes[-1].set_xlabel('Time [s]')

fig.tight_layout(rect=(0, 0, 1, 0.98))
fig.subplots_adjust(hspace=.0)