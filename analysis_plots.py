#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:04:43 2020

@author: misiak
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

from representation_functions import (
    temporal_plot,
    temporal_plot_heat_only,
    plot_chi2_vs_energy,
    histogram_adu,
    histogram_ev,
    ion_vs_ion,
    virtual_vs_virtual_ev,
    crosstalk_correction,
    trigger_cut_plot,
    plot_10kev,
    fid_cut_plot,
    band_cut_plots,
    ## nodecor
    nodecor_crosstalk_correction,
    charge_conservation,
)


def save_figure_dict(fig_dict, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for key, fig in fig_dict.items():
        fig.savefig( output_dir + '/{}.png'.format(key) )


def precalibration_plots(
        stream,
        title,
        df_analysis
    ):
    
    fig_dict = dict()
    
    ### temporal multi
    fig_temp = temporal_plot(title, df_analysis)
    fig_dict['temporal_multi'] = fig_temp
    
    ### temporal heat
    fig_temp_heat = temporal_plot_heat_only(title, df_analysis)
    fig_dict['temporal_heat'] = fig_temp_heat
    
    ### chi2 plot
    fig_chi2 = plot_chi2_vs_energy(title, df_analysis, stream)
    fig_dict['chi2_plot'] = fig_chi2

    ### histogramm ADU
    fig_hist_trig = histogram_adu(title, df_analysis, bins=10000)
    fig_dict['histogramm_ADU'] = fig_hist_trig

    ### crosstalk correction
    fig_cross = crosstalk_correction(title, df_analysis)
    fig_dict['crosstalk_correction'] = fig_cross

    # ## nodecor
    # ### crosstalk correction
    # fig_cross_nodecor = nodecor_crosstalk_correction(title, df_analysis)
    # fig_dict['nodecor_crosstalk_correction'] = fig_cross_nodecor

    return fig_dict


def nodecor_plots(
        title,
        df_analysis
    ):
    
    fig_dict = dict()
    
    ## nodecor
    ### crosstalk correction
    fig_cross_nodecor = nodecor_crosstalk_correction(title, df_analysis)
    fig_dict['nodecor_crosstalk_correction'] = fig_cross_nodecor

    ### plot charge conservation
    fig_charge = charge_conservation(title, df_analysis)
    fig_dict['charge_conservation'] = fig_charge
    
    return fig_dict


def simu_only_plots(
        title,
        df_analysis
    ):
    
    fig_dict = {
        'trigger_cut': trigger_cut_plot(title, df_analysis)
    }

    return fig_dict
    

def postcalibration_plots(
        title,
        df_analysis
    ):
    
    fig_dict = dict()
    
    ### histogramm ev
    fig_hist_trig_ev = histogram_ev(title, df_analysis, bins=10000)
    fig_dict['histogramm_ev'] = fig_hist_trig_ev
    
    ### ion vs ion
    fig_ion = ion_vs_ion(title, df_analysis)
    fig_dict['ion_vs_ion'] = fig_ion
    
    ### ion vs ion virtual
    fig_virtual = virtual_vs_virtual_ev(title, df_analysis)
    fig_dict['ion_vs_ion_virtual'] = fig_virtual
    
    ### 10kev plot
    fig_10kev = plot_10kev(title, df_analysis)
    fig_dict['plot_10kev'] = fig_10kev

    ### plot fig cut
    fig_fid = fid_cut_plot(title, df_analysis)
    fig_dict['fid_cut'] = fig_fid

    ### plot band cut
    fig_ecei, fig_quenching = band_cut_plots(title, df_analysis)
    fig_dict['band_cut_ecei'] = fig_ecei
    fig_dict['band_cut_quenching'] = fig_quenching
    
    # ### plot charge conservation
    # fig_charge = charge_conservation(title, df_analysis)
    # fig_dict['charge_conservation'] = fig_charge
    
    return fig_dict


def plotting_data_stream(
        stream,
        hdf5_path,
        output_dir,
        save_flag=False
    ):

    # saving all the figures
    save_dir = '/'.join([
        output_dir,
        stream, 
    ])
    
    df_analysis = pd.read_hdf(
        hdf5_path,
        key='df',
        where=(
            'stream = "{0}"'
        ).format(stream)
    )

    source = df_analysis['source'].unique()[0]

    title = '{0} {1} {2}'.format(stream, 'data', source).replace('_', ' ')

    fig_dict = {
        **precalibration_plots(stream, title, df_analysis),
        **postcalibration_plots(title, df_analysis),
        **nodecor_plots(title, df_analysis),
    }
    
    if save_flag:
        save_figure_dict(fig_dict, save_dir)
        
    return fig_dict


def plotting_noise_stream(
        stream,
        hdf5_path,
        output_dir,
        save_flag=False
    ):

    # saving all the figures
    save_dir = '/'.join([
        output_dir,
        stream,
        'noise'
    ])
    
    df_analysis = pd.read_hdf(
        hdf5_path,
        key='df',
        where=(
            'stream = "{0}"'
        ).format(stream)
    )

    source = df_analysis['source'].unique()[0]

    title = '{0} {1} {2}'.format(stream, 'noise', source).replace('_', ' ')

    fig_dict = {
        **precalibration_plots(stream, title, df_analysis),
        **postcalibration_plots(title, df_analysis),
    }
    
    if save_flag:
        save_figure_dict(fig_dict, save_dir)
        
    return fig_dict


def plotting_simu_stream(
        stream,
        simulation,
        hdf5_path,
        output_dir,
        save_flag=False
    ):

    # saving all the figures
    save_dir = '/'.join([
        output_dir,
        stream,
        simulation
    ])
    
    df_analysis = pd.read_hdf(
        hdf5_path,
        key='df',
        where=(
            'stream = "{0}"'
            '& simulation = "{1}"'
        ).format(stream, simulation)
    )

    source = df_analysis['source'].unique()[0]

    title = '{0} {1} {2}'.format(stream, simulation, source).replace('_', ' ')

    fig_dict = {
        **precalibration_plots(stream, title, df_analysis),
        **simu_only_plots(title, df_analysis),
        **postcalibration_plots(title, df_analysis),
        **nodecor_plots(title, df_analysis),
    }
    
    if save_flag:
        save_figure_dict(fig_dict, save_dir)
        
    return fig_dict


def plotting_data_source(
        source,
        hdf5_path,
        output_dir,
        save_flag=False
    ):

    # saving all the figures
    save_dir = '/'.join([
        output_dir,
        source, 
    ])
    
    df_analysis = pd.read_hdf(
        hdf5_path,
        key='df',
        where=(
            'source = "{0}"'
        ).format(source)
    )

    title = '{0} {1} {2}'.format('All streams', 'Data', source).replace('_', ' ')

    fig_dict = {
        **postcalibration_plots(title, df_analysis),
        **nodecor_plots(title, df_analysis),
    }
    
    if save_flag:
        save_figure_dict(fig_dict, save_dir)
        
    return fig_dict


def plotting_noise_source(
        source,
        hdf5_path,
        output_dir,
        save_flag=False
    ):

    # saving all the figures
    save_dir = '/'.join([
        output_dir,
        source,
        'noise'
    ])
    
    df_analysis = pd.read_hdf(
        hdf5_path,
        key='df',
        where=(
            'source = "{0}"'
        ).format(source)
    )

    title = '{0} {1} {2}'.format('All streams', 'Noise', source).replace('_', ' ')

    fig_dict = {
        **postcalibration_plots(title, df_analysis),
    }
    
    if save_flag:
        save_figure_dict(fig_dict, save_dir)
        
    return fig_dict


def plotting_simu_source(
        source,
        simulation,
        hdf5_path,
        output_dir,
        save_flag=False
    ):

    # saving all the figures
    save_dir = '/'.join([
        output_dir,
        source,
        simulation
    ])
    
    df_analysis = pd.read_hdf(
        hdf5_path,
        key='df',
        where=(
            'source = "{0}"'
            '& simulation = "{1}"'
        ).format(source, simulation)
    )

    source = df_analysis['source'].unique()[0]

    title = '{0} {1} {2}'.format("All streams", simulation, source).replace('_', ' ')

    fig_dict = {
        **simu_only_plots(title, df_analysis),
        **postcalibration_plots(title, df_analysis),
        **nodecor_plots(title, df_analysis),
    }
    
    if save_flag:
        save_figure_dict(fig_dict, save_dir)
        
    return fig_dict


if __name__ == "__main__":

    plt.close('all')
    from tqdm import tqdm
    
    analysis_dir = '/home/misiak/Analysis/neutron_background'
    output_plot_dir = '/'.join([analysis_dir, 'analysis_plots'])
    
    analysis_data_path = '/'.join([analysis_dir, 'data_analysis.h5'])
    analysis_noise_path = '/'.join([analysis_dir, 'noise_analysis.h5'])
    analysis_simu_path = '/'.join([analysis_dir, 'simu_analysis.h5'])
    
    stream_list = [
        'tg18l005',
        'tg27l000',
        'tg28l000',
        'tg17l007',
        'tg19l010',
        'tg20l000',
        'tg21l000'
    ]

    source_list = ['Calibration', 'Background']
    
    simulation_list = [
        'flat_ER',
        'flat_NR',
        'line_1keV',
        'line_10keV',
    ]

    # plotting_noise_stream(
    #     'tg28l000',
    #     analysis_noise_path,
    #     output_plot_dir,
    #     save_flag=False
    # )

    # A = plotting_noise_source(
    #     'Background',
    #     analysis_noise_path,
    #     output_plot_dir,
    #     save_flag=True
    # )

    A = plotting_data_source(
        'Background',
        analysis_data_path,
        output_plot_dir,
        save_flag=False
    )

    # for stream in tqdm(stream_list):
    #     plotting_data_stream(
    #         stream,
    #         analysis_data_path,
    #         output_plot_dir,
    #         save_flag=True
    #     )
    #     plt.close('all')
    # print('Stream Data done.')

    # for source in tqdm(source_list):
    #     plotting_data_source(
    #         source,
    #         analysis_data_path,
    #         output_plot_dir,
    #         save_flag=True
    #     )
    #     plt.close('all')
    # print('Source Data done.')

    # for stream in tqdm(stream_list):
    #     for simulation in simulation_list:
    #         plotting_simu_stream(
    #             stream,
    #             simulation,
    #             analysis_simu_path,
    #             output_plot_dir,
    #             save_flag=True
    #         )
    #         plt.close('all')
    # print('Stream Simulation done.')
 
    # for source in tqdm(source_list):
    #     for simulation in simulation_list:
    #         plotting_simu_source(
    #             source,
    #             simulation,
    #             analysis_simu_path,
    #             output_plot_dir,
    #             save_flag=True
    #         )
    #         plt.close('all')
    # print('Source Simulation done.')
        