#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 14:55:53 2020

@author: misiak
"""

import numpy as np
import pandas as pd
from tqdm import tqdm

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
            [7, 11.3],
            [13.8, 14.1]
        ],
        "heat_chi2_threshold": 700,
        "position_10kev_line_heat": 1293,
    },
    "tg28l000": {
    	"glitch_time_cut": [ 
            [7.4, 8.05],
            [10, np.inf]
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
    	"glitch_time_cut": [
            [8.3, 8.7]    
        ],
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

nodecor_crosstalk_correction_matrix = np.array([
        [1, -0.044, 0, 0],
        [-0.03, 1, 0, 0],
        [-0.025, 0.001, 1, -0.031],
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
        
        truth_array = pd.Series(data=True, index=df.index)
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
    
    truth_array = pd.Series(data=True, index=df.index)
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
    
    truth_array = pd.Series(data=True, index=df.index)
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


def nodecor_crosstalk_correction(df):
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


def nodecor_calibration_ion(stream, df):
    """ 
    Create new columns for the calibrated energy of the ionization channels.
    """         
    for suffix in ion_channel_labels:
        
        position_cname = 'position_10kev_line_ion{}'.format(suffix)
        energy_adu_cname ='energy_adu_corr_nodecor_ion{}'.format(suffix)
        energy_cname = 'energy_nodecor_ion{}'.format(suffix)
        
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


def nodecor_virtual_channels(df):
    """ 
    Create new columns for "virtual channels" which are combinations
    of the ionization channels:
        - energy_nodecor_ion_total: A+B+C+D
        - energy_nodecor_ion_collect: B+D
        - energy_nodecor_ion_guard: A+C
    """    

    df['energy_nodecor_ion_total'] = (
        df['energy_nodecor_ionA'] + df['energy_nodecor_ionB']
        + df['energy_nodecor_ionC'] + df['energy_nodecor_ionD']
    ) / 2
    
    df['energy_nodecor_ion_bulk'] = (
        df['energy_nodecor_ionB'] + df['energy_nodecor_ionD']
    ) / 2
    
    df['energy_nodecor_ion_guard'] = (
        df['energy_nodecor_ionA'] + df['energy_nodecor_ionC']
    ) / 2
    
    # XXX hard coded polarization
    df['energy_nodecor_ion_conservation'] = (
        - df['energy_nodecor_ionA'] - df['energy_nodecor_ionB']
        + df['energy_nodecor_ionC'] + df['energy_nodecor_ionD']
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

def energy_heat_from_er_and_quenching(er, Q, V): 
    return er * (1 + Q*V/3) / (1 + V/3)

def energy_ion_from_er_and_quenching(er, Q, V=None): 
    return er * Q


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


def guard_threshold_for_bulk_cut(energy_ion):
    # std_noise_blob_fid = 1
    # std_10kev_fid = 1.5
    std_noise_blob_fid = 0.29
    std_10kev_fid = 0.5
    alpha_fid = (std_10kev_fid**2 - std_noise_blob_fid**2)**0.5 / 10.37    
    return (std_noise_blob_fid**2 + (alpha_fid*energy_ion)**2)**0.5


def bulk_threshold_for_guard_cut(energy_ion):
    # std_noise_blob_fid = 1
    # std_10kev_fid = 1.5
    std_noise_blob_fid = 0.29
    std_10kev_fid = 0.5
    alpha_fid = (std_10kev_fid**2 - std_noise_blob_fid**2)**0.5 / 10.37    
    return (std_noise_blob_fid**2 + (alpha_fid*energy_ion)**2)**0.5


def fid_cuts(df):
    """ 
    Apply the so-called FID cuts, which is a way to discriminate events
    happening in the bulk (or in the guard) region from the others.
    Create new columns with the truth array for the bulk and guard events.
    """    
    energy_ion_collect = df['energy_ion_bulk']
    energy_ion_guard = df['energy_ion_guard']
    
    nsigma = 2
    
    collect_energy_threshold = guard_threshold_for_bulk_cut(energy_ion_collect) * nsigma

    guard_energy_threshold = bulk_threshold_for_guard_cut(energy_ion_guard) * nsigma
    
    df['bulk_cut'] = (
        ( abs(df['energy_ionA']) < collect_energy_threshold )
        & ( abs(df['energy_ionC']) < collect_energy_threshold )
    )

    df['guard_cut'] = (
        ( abs(df['energy_ionB']) < guard_energy_threshold )
        & ( abs(df['energy_ionD']) < guard_energy_threshold )
    )
    
    return None


def charge_conservation_threshold(energy_heat):
    nsigma = 2
    return (0.4120323 + 0.00165406 * energy_heat) * nsigma

def charge_conservation_cut(df):
    """ 
    Create a new column with the truth array of the 
    charge conservation cut.
    """    
    ion_conservation = df['energy_nodecor_ion_conservation']
    energy_heat = df['energy_heat']
    
    threshold = charge_conservation_threshold(energy_heat)
        
    df['charge_conservation_cut'] = abs(ion_conservation) < threshold
        
    return None


# def nodecor_fid_cuts(df):
#     """ 
#     Apply the so-called FID cuts, which is a way to discriminate events
#     happening in the bulk (or in the guard) region from the others.
#     Create new columns with the truth array for the bulk and guard events.
#     """    
#     energy_ion_collect = df['energy_nodecor_ion_bulk']
#     energy_ion_guard = df['energy_nodecor_ion_guard']
    
#     collect_energy_threshold = guard_threshold_for_bulk_cut(energy_ion_collect)

#     guard_energy_threshold = bulk_threshold_for_guard_cut(energy_ion_guard)
    
#     df['nodecor_bulk_cut'] = (
#         ( abs(df['energy_nodecor_ionA']) < collect_energy_threshold )
#         & ( abs(df['energy_nodecor_ionC']) < collect_energy_threshold )
#     )

#     df['nodecor_guard_cut'] = (
#         ( abs(df['energy_nodecor_ionB']) < guard_energy_threshold )
#         & ( abs(df['energy_nodecor_ionD']) < guard_energy_threshold )
#     )
    
#     return None


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


def trigger_cut(
        df,
        nearest_data_trigger_allowed=5,
        furthest_simu_trigger_allowed=5
    ):
    
    df['trigger_cut'] = (
        ( abs(df['t_nearest_data_trigger']) > nearest_data_trigger_allowed )
        & ( abs(df['t_input_simu_trigger']) < furthest_simu_trigger_allowed )
    )


def std_energy_ion(ec):
    """
    standard deviation on the ionization energy in function of the heat energy.
    """
    # std_nb = 0.254
    # std_nb = 0.200
    std_nb = 0.22
    #std_1kev = 0.272
    # std_10kev = 0.318
    std_10kev = 0.45
    alpha = (std_10kev**2 - std_nb**2)**0.5 / 10.37
    return ( std_nb**2 + (alpha*ec)**2 )**0.5


def band_cut(df):

    energy_heat = df['energy_heat']
    energy_ion = df['energy_ion_bulk']
    nsigma = 2
    ei_err = nsigma * std_energy_ion(energy_heat)

    gamma_cut = ( abs(energy_ion - energy_heat) < ei_err )

    er_array = np.linspace(0, energy_heat.max(), int(1e4))
    dv=2
    ec_array = er_array * (1 + lindhard(er_array)*dv/3) / (1 + dv/3)
    ei_array = er_array * lindhard(er_array)
    energy_ion_lindhard = np.interp(energy_heat, ec_array, ei_array)
    
    neutron_cut = ( abs(energy_ion - energy_ion_lindhard) < ei_err )

    HO_cut = ( abs(energy_ion) < ei_err )

    df['gamma_cut'] = gamma_cut
    df['neutron_cut'] = neutron_cut
    df['HO_cut'] = HO_cut

    return None


def analysis_data(stream, df_fine):
    """
    Use with experimental data.
    Create an analysed dataframe from the given fine dataframe.
    """
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
    band_cut(df_analysis)

    #nodecor
    nodecor_crosstalk_correction(df_analysis)
    nodecor_calibration_ion(stream, df_analysis)
    
    nodecor_virtual_channels(df_analysis)
    charge_conservation_cut(df_analysis)
    # nodecor_fid_cuts(df_analysis)    

    return df_analysis


def analysis_noise(stream, df_fine):
    """
    Use with experimental data.
    Create an analysed dataframe from the given fine dataframe.
    """
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
    band_cut(df_analysis)  

    return df_analysis


def analysis_simulation(stream, df_fine):
    """
    Use with pulse simulation.
    Create an analysed dataframe from the given fine dataframe.
    """
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
    trigger_cut(df_analysis)
    energy_cut(df_analysis)
    band_cut(df_analysis)
    
    physical_quantities(df_analysis)

    #nodecor
    nodecor_crosstalk_correction(df_analysis)
    nodecor_calibration_ion(stream, df_analysis)
    
    nodecor_virtual_channels(df_analysis)
    charge_conservation_cut(df_analysis)
    # nodecor_fid_cuts(df_analysis)    

    
    return df_analysis


def hdf5_analysis(fine_hdf5_path, output_hdf5_path, analysis_function):

    stream_list = pd.read_hdf(
        fine_hdf5_path,
        key='df',
        columns=['stream',]
    )['stream'].unique()
    
    # initializing the HDFstore (overwriting, be careful !)
    pd.DataFrame().to_hdf(
        output_hdf5_path,
        key='df', mode='w', format='table'
    )
    
    for stream in tqdm(stream_list):

        df_fine = pd.read_hdf(
            fine_hdf5_path,
            key='df',
            where='stream = "{}"'.format(stream)
        )
        
        df_analysis = analysis_function(stream, df_fine)

        df_analysis.to_hdf(
            output_hdf5_path,
            key='df',
            mode='a',
            format='table',
            append=True,
            min_itemsize=11,
            data_columns=True
        ) 

    return None
    

if __name__ == "__main__":

    analysis_dir = '/home/misiak/Analysis/neutron_background'
    
    fine_data_path = '/'.join([analysis_dir, 'data_fine.h5'])
    output_data_path = '/'.join([analysis_dir, 'data_analysis.h5'])   
    
    hdf5_analysis(
        fine_data_path,
        output_data_path,
        analysis_data
    )

    fine_noise_path = '/'.join([analysis_dir, 'noise_fine.h5'])
    output_noise_path = '/'.join([analysis_dir, 'noise_analysis.h5'])   
    
    hdf5_analysis(
        fine_noise_path,
        output_noise_path,
        analysis_noise
    )
    
    fine_simu_path = '/'.join([analysis_dir, 'simu_fine.h5'])
    output_simu_path = '/'.join([analysis_dir, 'simu_analysis.h5'])
    hdf5_analysis(
        fine_simu_path,
        output_simu_path,
        analysis_simulation
    )
    