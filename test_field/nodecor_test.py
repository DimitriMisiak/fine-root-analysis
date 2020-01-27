#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 11:52:40 2020

@author: misiak
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from representation_functions import nodecor_crosstalk_correction, ax_hist
from data_analysis import (
    glitch_time_cut,
    heat_chi2_cut,
    ion_chi2_cut,
    offset_ion_cut,
    quality_cut
)
    

fine_data_path = '/home/misiak/Analysis/neutron_background/data_fine.h5'

stream = 'tg20l000'
DF = pd.read_hdf(
    fine_data_path,
    where='stream = "{}"'.format(stream))

glitch_time_cut(stream, DF)
heat_chi2_cut(stream, DF)
ion_chi2_cut(stream, DF)
offset_ion_cut(DF)
quality_cut(DF)

#%%
plt.close('all')

nodecor_crosstalk_correction_matrix = np.array([
        [1, -0.044, 0, 0],
        [-0.03, 1, 0, 0],
        [-0.025, 0.001, 1, -0.031],
        [0, 0, -0.03, 1]
])

def funk(df):
    """ 
    Create new columns for the cross-talk corrected ionization channels.
    """        
    
    ion_energy = df[[
        'energy_adu_nodecor_ionA',
        'energy_adu_nodecor_ionB',
        'energy_adu_nodecor_ionC',
        'energy_adu_nodecor_ionD'
    ]]
    
    
    corr_ion_energy = np.dot(nodecor_crosstalk_correction_matrix, ion_energy.T)
    energy_corr_cname_list = [
        'energy_adu_corr_nodecor_ionA',
        'energy_adu_corr_nodecor_ionB',
        'energy_adu_corr_nodecor_ionC',
        'energy_adu_corr_nodecor_ionD'
    ]
    for i, col in enumerate(energy_corr_cname_list):
        df[col] = corr_ion_energy[i]
    
    return None

funk(DF)

nodecor_crosstalk_correction('hell', DF)


def histogram_adu(title, df, bins=1000):
    
    ax_tuples = ((0, 1), (1, 0), (1, 1), (2, 0), (2, 1))
 
    channel_suffix = [
        'heat',
        'ionA',
        'ionB',
        'ionC',
        'ionD',
    ]    
    
    ax_titles =[
        'Heat',
        'Ion A',
        'Ion B',
        'Ion C',
        'Ion D',
    ]
    
    quality_cut = df['quality_cut']
    
    num = '{} : Quality Cut Histogram'.format(title)

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(11.69, 8.27),
                             num=num)
    
    for suffix, tupl, label in zip(channel_suffix, ax_tuples, ax_titles):
        
        if suffix == 'heat':
            xdata = df['energy_adu_{}'.format(suffix)]
        else:
            xdata = df['energy_adu_corr_nodecor_{}'.format(suffix)]
        
        ax = axes[tupl]
        xdata_qual = xdata[quality_cut]
        

        # try:
        #     bin_edges = custom_bin_edges(xdata_qual, 
        #                                  getattr(noise.sigma0, label))
        
        bin_edges = np.histogram_bin_edges(xdata[quality_cut], bins=bins)

        ax_hist(ax, bin_edges, xdata,
                'All events', color='coral')
        ax_hist(ax, bin_edges, xdata_qual,
                'Quality events', color='slateblue')[0]
        


        xdata_fid = xdata[quality_cut]
        ax_hist(ax, bin_edges, xdata_fid,
                'Fiducial events', color='limegreen')[0]        
    

        ax.legend(loc=2)
        ax.set_title(label.replace('_', ' '))
        
    
    fig.text(0.5, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))

    # resize the plots
    fig.get_axes()[0].set_xlim(-200, 2000)
    for i, ax in enumerate(fig.get_axes()[:5]):
        if i==0:
            ax.set_xlim(-200, 2000)
        else:
            ax.set_xlim(-70, 70)    
 
    fig.delaxes(axes[0,0])    
    fig.tight_layout()
        
    return fig

histogram_adu('hell yeah', DF)
