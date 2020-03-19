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

fig, axes = plt.subplots(ncols=2)

D = df_analysis[quality_cut & C]

xB = D['energy_adu_corr_ionB']
yA = D['energy_adu_corr_ionA']
line, = axes[0].plot(
    xB,
    yA,
    ls='none',
    marker='.'
)

axes[0].set_xlabel('Amplitude B')
axes[0].set_ylabel('Amplitude A')
axes[0].grid()


linear = lambda x,*p: p[0]*x + p[1]
from scipy.optimize import curve_fit


def funAB(indexes):
    x_data = D['energy_adu_ionB'].array[indexes]
    y_data = D['energy_adu_ionA'].array[indexes]

    popt, pcov = curve_fit(
        linear,
        x_data,
        y_data,
        p0=[1,1]
    )

    ion_array = np.linspace(
        x_data.min(),
        x_data.max(),
        int(1e3)
    )       
    
    axes[0].plot(
        ion_array,
        linear(ion_array, *popt),
        color='k'
    )
    
    # axes[0].legend(title=popt)
    print(popt)
    return popt
    
# ds = Data_Selector(axes[0], line, funAB)

from scipy.odr import ODR, Model, Data, RealData

# def func(beta, x):
#     y = beta[0]+beta[1]*x
#     return y

def func(beta, x):
    y = (-beta[1]/beta[0])*x+beta[1]
    return y


### Calibration A & B
data = RealData(xB, yA, 5, 5)
model = Model(func)
   
odr = ODR(data, model, [1,1])

odr.set_job(fit_type=0)
output = odr.run()

xn = np.linspace(
    xB.min(),
    xB.max(),
    int(1e3)
)       
yn = func(output.beta, xn)

axes[0].plot(xn,yn,'r-',label='odr')


### Calibration C & D

xD = D['energy_adu_corr_ionD']
yC = D['energy_adu_corr_ionC']

axes[1].plot(
    xD,
    yC,
    ls='none',
    marker='.'
)

axes[1].set_xlabel('Amplitude D')
axes[1].set_ylabel('Amplitude C')
axes[1].grid()

data = RealData(xD, yC, 5, 5)
model = Model(func)
   
odr = ODR(data, model, [1,1])

odr.set_job(fit_type=0)
output = odr.run()

xn = np.linspace(
    xD.min(),
    xD.max(),
    int(1e3)
)       
yn = func(output.beta, xn)

axes[1].plot(xn,yn,'r-',label='odr')


from scipy.stats import binned_statistic

lambda_median = lambda x: np.quantile(x, 0.5)
lambda_inf = lambda x: np.quantile(x, 0.16)
lambda_sup = lambda x: np.quantile(x, 0.84)
median, bin_edges, _ = binned_statistic(xD, yC, lambda_median)
inf, bin_edges, _ = binned_statistic(xD, yC, lambda_sup)
sup, bin_edges, _ = binned_statistic(xD, yC, lambda_inf)

bin_width = (bin_edges[1] - bin_edges[0])
bin_array = bin_edges[:-1] + bin_width/2

axes[1].errorbar(
    bin_array,
    median,
    yerr=[median-inf, sup-median],
    xerr=bin_width/2,
    marker='o',
    color='k'
)

# fig, ax = plt.subplots()

# bin_edges = np.linspace(1000, 1400, int(1e2))

# ax0 = ax_hist(ax, bin_edges, A[B], 'quality')[0]

# import scipy.stats as st
# popt = st.norm.fit(A[B])

# pdf = st.norm(*popt).pdf(bin_edges)
# cdf = st.norm(*popt).cdf(bin_edges)

# # ax.plot(bin_edges, pdf, color='k')
# ax0.plot(bin_edges, cdf, color='k')
    
    