#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:55:53 2020

@author: misiak
"""

import numpy as np
import pandas as pd


analysis_parameters = {
    "tg18l005": {
    	"glitch_time_cut": [ 
            [7.4, 7.6]
        ],
        "heat_chi2_threshold": 400,
        "position_10kev_line_heat": 1177,
    },
    "tg27l000": {
    	"glitch_time_cut": [ 
            [7, 11.3]
        ],
        "heat_chi2_threshold": 700,
        "position_10kev_line_heat": 1293,
    },
    "tg28l000": {
    	"glitch_time_cut": [ 
            [7.4, 8.05]
        ],
        "heat_chi2_threshold": 700,
        "position_10kev_line_heat": 1291,
    },
    "tg17l007": {
    	"glitch_time_cut": [],
        "heat_chi2_threshold": 400,
        "position_10kev_line_heat": 1165,
    },
    "tg19l010": {
    	"glitch_time_cut": [],
        "heat_chi2_threshold": 400,
        "position_10kev_line_heat": 1199,
    },
    "tg20l000": {
    	"glitch_time_cut": [],
        "heat_chi2_threshold": 400,
        "position_10kev_line_heat": 1197,
    },
    "tg21l000": {
    	"glitch_time_cut": [],
        "heat_chi2_threshold": 400,
        "position_10kev_line_heat": 1235,
    },
    "ion_chi2_threshold": 300,
    "offset_ion_threshold": 14000,
    "position_10kev_line_ionA": -56,
    "position_10kev_line_ionB": -54,
    "position_10kev_line_ionC": 57,
    "position_10kev_line_ionD": 56,
}

ion_channel_labels = ('A', 'B', 'C', 'D')

crosstalk_correction_matrix = np.array([
        [1, -0.052, 0, 0],
        [-0.03, 1, 0, 0],
        [-0.012, 0.001, 1, -0.025],
        [0, 0, -0.03, 1]
])


def glitch_time_cut(stream, df):
    """ 
    Create a new column with the "glitch time cut" truth array.
    """
    glitch_interval_list = analysis_parameters[stream]['glitch_time_cut']
    
    if glitch_interval_list is None:
        df['glitch_time_cut'] = True
        
    else:
        
        timestamp = df['timestamp']
        
        truth_array = pd.Series(data=True, index=df_fine.index)
        for inf, sup in glitch_interval_list:
            
            interval_truth = (timestamp < inf) | (timestamp > sup)
            truth_array = truth_array & interval_truth
            
        df['glitch_time_cut'] = truth_array
    
    return None
    

def heat_chi2_threshold_function(x2, e0):
    return x2 * ( 1 + (e0/2e3)**2 )

def heat_chi2_cut(stream, df):
    """ 
    Create a new column with the truth array of the 
    chi2 cut on the heat channel.
    """
    heat_chi2_threshold = analysis_parameters[stream]['heat_chi2_threshold']
    chi2_heat = df['chi2_heat']
    energy_adu_heat = df['energy_adu_heat']
    
    # chi2_threshold = heat_chi2_threshold * ( 1 + (energy_adu_heat/2e3)**2 )
    chi2_threshold = heat_chi2_threshold_function(
        heat_chi2_threshold,
        energy_adu_heat
    )
        
    df['chi2_heat_cut'] = (chi2_heat < chi2_threshold)
    
    return None

def ion_chi2_threshold_function(x2, e0):
    return x2 * ( 1 + (e0/3e2)**2.2 )

def ion_chi2_cut(stream, df):
    """ 
    Create a new column with the truth array of the 
    chi2 cut on the ion channel.
    """
    ion_chi2_threshold = analysis_parameters['ion_chi2_threshold']
    
    truth_array = pd.Series(data=True, index=df_fine.index)
    for suffix in ion_channel_labels:
        chi2_column_name = 'chi2_ion{}'.format(suffix)
        chi2_ion = df[chi2_column_name]
        
        energy_column_name = 'energy_adu_ion{}'.format(suffix)
        energy_adu_ion = abs( df[energy_column_name] )
        
        # chi2_threshold = ion_chi2_threshold * ( 1 + (energy_adu_ion/3e2)**2.2 )
        chi2_threshold = ion_chi2_threshold_function(
            ion_chi2_threshold,
            energy_adu_ion
        )
        
        truth_array = truth_array & (chi2_ion < chi2_threshold)
        
    df['chi2_ion_cut'] = truth_array
    
    return None


def offset_ion_cut(df):
    """ 
    Create a new column with the truth array of the 
    chi2 cuts on the ion channels.
    """    

    offset_threshold = analysis_parameters['offset_ion_threshold']
    
    truth_array = pd.Series(data=True, index=df_fine.index)
    for suffix in ion_channel_labels:
        offset_column_name = 'offset_ion{}'.format(suffix)
        offset = abs( df[offset_column_name] )
        truth_array = truth_array & (offset < offset_threshold)
        
    df['offset_ion_cut'] = truth_array
    
    return None

    
def quality_cut(df):
    """ 
    Create a new column with the truth array of the 
    quality cut.
    """    

    quality_cut_components_columns = [
        'glitch_time_cut',
        'offset_ion_cut',
        'chi2_heat_cut',
        'chi2_ion_cut',
    ]
    
    truth_array = pd.Series(data=True, index=df.index)
    for col in quality_cut_components_columns:
        truth_array = truth_array & df[col]
        
    df['quality_cut'] = truth_array

    return None



def crosstalk_correction(df):
    """ 
    Create new columns for the cross-talk corrected ionization channels.
    """        
    
    ion_energy = df[[
        'energy_adu_ionA',
        'energy_adu_ionB',
        'energy_adu_ionC',
        'energy_adu_ionD'
    ]]
    
    
    corr_ion_energy = np.dot(crosstalk_correction_matrix, ion_energy.T)
    energy_corr_cname_list = [
        'energy_adu_corr_ionA',
        'energy_adu_corr_ionB',
        'energy_adu_corr_ionC',
        'energy_adu_corr_ionD'
    ]
    for i, col in enumerate(energy_corr_cname_list):
        df[col] = corr_ion_energy[i]
    
    return None


def calibration_heat(stream, df):
    """ 
    Create a new column for the calibrated energy of heat channel.
    """       
    calibration_factor = (
        10.37 / analysis_parameters[stream]['position_10kev_line_heat']
    )
    
    df['energy_heat'] = (
        df['energy_adu_heat'] * calibration_factor
    )
    
    return None


def calibration_ion(stream, df):
    """ 
    Create new columns for the calibrated energy of the ionization channels.
    """         
    for suffix in ion_channel_labels:
        
        position_cname = 'position_10kev_line_ion{}'.format(suffix)
        energy_adu_cname ='energy_adu_corr_ion{}'.format(suffix)
        energy_cname = 'energy_ion{}'.format(suffix)
        
        calibration_factor = (
            10.37 / analysis_parameters[position_cname]
        )
        
        df[energy_cname] = (
            df[energy_adu_cname] * calibration_factor
        )
        
    return None


def virtual_channels(df):
    """ 
    Create new columns for "virtual channels" which are combinations
    of the ionization channels:
        - energy_ion_total: A+B+C+D
        - energy_ion_collect: B+D
        - energy_ion_guard: A+C
    """    

    df['energy_ion_total'] = (
        df['energy_ionA'] + df['energy_ionB']
        + df['energy_ionC'] + df['energy_ionD']
    ) / 2
    
    df['energy_ion_bulk'] = (
        df['energy_ionB'] + df['energy_ionD']
    ) / 2
    
    df['energy_ion_guard'] = (
        df['energy_ionA'] + df['energy_ionC']
    ) / 2
    
    # XXX hard coded polarization
    df['energy_ion_conservation'] = (
        - df['energy_ionA'] - df['energy_ionB']
        + df['energy_ionC'] + df['energy_ionD']
    ) / 2

    return None


def energy_recoil(ec, ei, V):
#    coeff = 1.6e-19 * V / 3
    coeff = V / 3
    return ec*(1+coeff) - ei*coeff


def quenching(ec, ei, V):
    er = energy_recoil(ec, ei, V)
    return ei/er


def lindhard(er):
    
    A = 72.63
    Z = 32
    
    k = 0.133 * Z**(2./3) * A**(-1./2)
    epsilon = 11.5 * er * Z**(-7./3)
    g = 3 * epsilon**0.15 + 0.7 * epsilon**0.6 + epsilon  
    
    Q = k*g/(1+k*g)
    
    return Q


def physical_quantities(df, voltage=2):
    """ 
    Create new columns for the recoil energy and quenching.
    """    
    for suffix in ('_total', '_bulk', '_guard'):
        
        energy_ion_cname = 'energy_ion{}'.format(suffix)
        recoil_energy_cname = 'recoil_energy{}'.format(suffix)
        quenching_cname = 'quenching{}'.format(suffix)
        
        df[recoil_energy_cname] = energy_recoil(
            df['energy_heat'],
            df[energy_ion_cname],
            voltage
        )
    
        df[quenching_cname] = quenching(
            df['energy_heat'],
            df[energy_ion_cname],
            voltage
        )

    return None

    
def fid_cuts(df):
    """ 
    Apply the so-called FID cuts, which is a way to discriminate events
    happening in the bulk (or in the guard) region from the others.
    Create new columns with the truth array for the bulk and guard events.
    """    
    std_noise_blob_fid = 1
    
    std_10kev_fid = 1.5
    
    alpha_fid = (std_10kev_fid**2 - std_noise_blob_fid**2)**0.5 / 10.37
    
    energy_ion_collect = df['energy_ion_bulk']
    energy_ion_guard = df['energy_ion_guard']
    
    collect_energy_threshold = (
        (std_noise_blob_fid**2 + (alpha_fid*energy_ion_collect)**2)**0.5
    )

    guard_energy_threshold = (
        (std_noise_blob_fid**2 + (alpha_fid*energy_ion_guard)**2)**0.5
    )
    
    df['bulk_cut'] = (
        ( abs(df['energy_ionA']) < collect_energy_threshold )
        & ( abs(df['energy_ionC']) < collect_energy_threshold )
    )

    df['guard_cut'] = (
        ( abs(df['energy_ionB']) < guard_energy_threshold )
        & ( abs(df['energy_ionD']) < guard_energy_threshold )
    )
    
    return None


def energy_cut(df, energy_bounds=[0.025, 50]):
    """
    Extra cut to select a specific range in energy.
    Quite handy to drop event with negative energy (not physical) and event 
    with a high energy (non-linearity becomes problematic).
    """
    inf, sup = energy_bounds
    
    df['energy_cut'] = (
        ( df['energy_heat'] > inf )
        & ( df['energy_heat'] < sup )
    )
    
    return None



if __name__ == "__main__":

    analysis_dir = '/home/misiak/Analysis/neutron_background'
    fine_data_path = '/'.join([analysis_dir, 'data_fine.h5'])
    analysis_data_path = '/'.join([analysis_dir, 'data_analysis.h5'])


    stream_list = pd.read_hdf(
        fine_data_path,
        key='df',
        columns=['stream',]
    )['stream'].unique()
    
    # initializing the HDFstore (overwriting, be careful !)
    pd.DataFrame().to_hdf(
        analysis_data_path,
        key='df', mode='w', format='table'
    )
    
    from tqdm import tqdm
    for stream in tqdm(stream_list):

        df_fine = pd.read_hdf(
            fine_data_path,
            key='df',
            where='stream = "{}"'.format(stream)
        )
        
        df_analysis = pd.DataFrame()
        for col in df_fine.columns:
            df_analysis[col] = df_fine[col]
            
        glitch_time_cut(stream, df_analysis)
        heat_chi2_cut(stream, df_analysis)
        ion_chi2_cut(stream, df_analysis)
        offset_ion_cut(df_analysis)
        quality_cut(df_analysis)
        
        crosstalk_correction(df_analysis)
        calibration_heat(stream, df_analysis)
        calibration_ion(stream, df_analysis)
        
        virtual_channels(df_analysis)
        fid_cuts(df_analysis)
        energy_cut(df_analysis)
        
        physical_quantities(df_analysis)

        df_analysis.to_hdf(
            analysis_data_path,
            key='df',
            mode='a',
            format='table',
            append=True,
            min_itemsize=11,
            data_columns=True
        ) 
