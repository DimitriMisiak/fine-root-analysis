#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""

import matplotlib.pyplot as plt
from spec_classes import Analysis_red
from representation import (
        temporal_plot, plot_chi2_vs_energy,
        histogram_adu, histogram_ev, ion_vs_ion,
        virtual_vs_virtual_ev, optimization_info
)

run_dir = '/home/misiak/Data/data_run57'
stream = 'tg27l000'
detector = 'RED80'


# first command
plt.close('all')
plt.rcParams['text.usetex']=True


ana = Analysis_red(
        stream,
        detector=detector,
        run_dir=run_dir,
        chan_valid=(0, 2, 3, 4, 5),
        chan_signal=(0,)
)   

fig_temp = temporal_plot(ana)
fig_chi2_trig, fig_chi2_noise = plot_chi2_vs_energy(ana)
fig_hist_trig, fig_hist_noise = histogram_adu(ana)
fig_hist_trig_ev, fig_hist_noise_ev = histogram_ev(ana)
fig_ion = ion_vs_ion(ana)
fig_virtual = virtual_vs_virtual_ev(ana)

plt.show()
