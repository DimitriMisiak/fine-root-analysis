#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:04:43 2020

@author: misiak
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from plot_addon import LegendTitle, custom_autoscale, ax_hist


def temporal_plot(stream, df):
    """
    Monitoring plots.
    Several key quantities in function of time.

    Parameters
    ----------
    stream : str
        Stream name.
    df : pandas.DataFrame
        DF containing the analysed data from the "data_analysis.py" script.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    cut = df['quality_cut']
    
    source = df['source'].unique()[0]
    
    time = df['timestamp']
    energy_adu_heat = df['energy_adu_heat']
    chi2_heat = df['chi2_heat']
    offset_heat = df['offset_heat']
    
    energy_adu_ion = {
        'ionA': df['energy_adu_ionA'],
        'ionB': df['energy_adu_ionB'],
        'ionC': df['energy_adu_ionC'],
        'ionD': df['energy_adu_ionD'],        
    }
    
    offset_ion = {
        'ionA': df['offset_ionA'],
        'ionB': df['offset_ionB'],
        'ionC': df['offset_ionC'],
        'ionD': df['offset_ionD'],        
    }
    
    slope_ion = {
        'ionA': df['slope_ionA'],
        'ionB': df['slope_ionB'],
        'ionC': df['slope_ionC'],
        'ionD': df['slope_ionD'],        
    }    
    
    # Init figure
    num = '{1} ({0}): Monitoring'.format(source, stream)
    fig, axes = plt.subplots(nrows=6, ncols=1, figsize=(12, 10),
                             sharex=True, num=num)
    
    # heat trig vs time
    ax = axes[0]
    ax.set_ylabel('Energy Heat [ADU]')
    ax.set_yscale('symlog')
    
    ax.plot(
            time[cut], energy_adu_heat[cut],
            label='heat a', zorder=10,
            ls='none', marker='2', mew=0.8,
    )
    ax.autoscale(False)
    ax.plot(
            time, energy_adu_heat,
            label='All events',
            ls='none', marker=',', color='silver',
    )
    
    # ion trig vs time
    ax = axes[1]
    ax.set_ylabel('Energy Ion [ADU]')
    ax.set_yscale('symlog')
        
    for i, items in enumerate(energy_adu_ion.items()):
        lab, series = items
        ax.plot(
                time[cut], series[cut],
                label=lab, zorder=10,
                ls='none', marker=str(i+1), mew=0.8
        )        
        
        
    ax.autoscale(False)
    for lab, series in energy_adu_ion.items():
        ax.plot(
                time, series,
                label='All events',
                ls='none', marker=',', color='silver',
        )
    
    # # heat offset vs time
    ax = axes[2]
    ax.set_ylabel('Offset Heat [ADU]')    
    ax.plot(
            time[cut], offset_heat[cut],
            label='heat a', zorder=10,
            ls='none', marker='2', mew=0.8,
    )
    ax.autoscale(False)
    ax.plot(
            time, offset_heat,
            label='All events',
            ls='none', marker=',', color='silver',
    )
    
    # ion offset vs time
    ax = axes[3]
    ax.set_ylabel('Offset Ion [ADU]')
    
    for i, items in enumerate(offset_ion.items()):
        lab, series = items
        ax.plot(
                time[cut], series[cut],
                label=lab, zorder=10,
                ls='none', marker=str(i+1), mew=0.8
        )        
        
    ax.autoscale(False)
    for lab, series in offset_ion.items():
        ax.plot(
                time, series,
                label='All events',
                ls='none', marker=',', color='silver',
        )
    
    # ion slope vs time
    ax = axes[4]
    ax.set_ylabel('Slope Ion [ADU/s]')
    for i, items in enumerate(slope_ion.items()):
        lab, series = items
        ax.plot(
                time[cut], series[cut],
                label=lab, zorder=10,
                ls='none', marker=str(i+1), mew=0.8
        )        
        
    ax.autoscale(False)
    for lab, series in slope_ion.items():
        ax.plot(
                time, series,
                label='All events',
                ls='none', marker=',', color='silver',
        )
    
    # chi2 vs time
    ax = axes[5]
    ax.set_ylabel('$\chi^2$')
    ax.set_yscale('log')
    label = 'chi2 heat A'
    ax.plot(
            time[cut], chi2_heat[cut],
            label='heat a', zorder=10,
            ls='none', marker='2', mew=0.8,
    )
    ax.autoscale(False)
    ax.plot(
            time, chi2_heat,
            label='All events',
            ls='none', marker=',', color='silver',
    )
    
    
    # formatting the axes
    for ax in axes:
        ax.grid(True, alpha=0.3)
        
        # custom legend
        handles = ['Quality events:',]
        labels = ['',]
        for line in ax.get_lines():
            label = line.get_label()
            if label == 'All events':
                if label != labels[0]:
                    handles.insert(0, line)
                    labels.insert(0, label)
            else:
                handles.append(line)
                labels.append(label)
        
        # handler_map + LegendTitle allow for subtitle in legend
        ax.legend(
                handles, labels, loc=2, framealpha=1,
                bbox_to_anchor=(1.05, 1), borderaxespad=0.,
                handler_map={str: LegendTitle()}
        )
    
        if ax is not axes[-1]:
            # removing the first tick label
            yticks = ax.yaxis.get_major_ticks()
            yticks[0].label1.set_visible(False)
    
        if ax is axes[-1]:
            ax.set_xlabel('Time [hours]')
    
    
    fig.text(0.5, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))
    
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.subplots_adjust(hspace=.0)
    
    return fig


def temporal_plot_heat_only(stream, df):

    cut = df['quality_cut']
    
    source = df['source'].unique()[0]
    
    time = df['timestamp']
    energy_adu_heat = df['energy_adu_heat']
    
    # Init figure
    num = '{1} ({0}): Monitoring heat channel'.format(source, stream)
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(12, 5),
                             sharex=True, num=num)
    
    axes = [axes,]
    
    # heat trig vs time
    ax = axes[0]
    ax.set_ylabel('Energy Heat [ADU]')
    ax.set_yscale('symlog')
    
    ax.plot(
            time[cut], energy_adu_heat[cut],
            label='heat a', zorder=10,
            ls='none', marker='2', mew=0.8,
    )
    ax.autoscale(False)
    ax.plot(
            time, energy_adu_heat,
            label='All events',
            ls='none', marker=',', color='silver',
    )
    
    # formatting the axes
    for ax in axes:
        ax.grid(True, alpha=0.3)
        
        # custom legend
        handles = ['Quality events:',]
        labels = ['',]
        for line in ax.get_lines():
            label = line.get_label()
            if label == 'All events':
                if label != labels[0]:
                    handles.insert(0, line)
                    labels.insert(0, label)
            else:
                handles.append(line)
                labels.append(label)
        
        # handler_map + LegendTitle allow for subtitle in legend
        ax.legend(
                handles, labels, loc=2, framealpha=1,
                bbox_to_anchor=(1.05, 1), borderaxespad=0.,
                handler_map={str: LegendTitle()}
        )
    
        if ax is not axes[-1]:
            # removing the first tick label
            yticks = ax.yaxis.get_major_ticks()
            yticks[0].label1.set_visible(False)
    
        if ax is axes[-1]:
            ax.set_xlabel('Time [hours]')
    
    
    fig.text(0.5, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))
    
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.subplots_adjust(hspace=.0)
    
    return fig


def plot_chi2_vs_energy(stream, df):

    channel_suffix = [
        'heat',
        'ionA',
        'ionB',
        'ionC',
        'ionD',
    ]
    
    # chi2 vs Energy plot
    ax_tuples = ((1, 0), (0, 1), (0, 2), (1, 1), (1, 2))
    
    ax_titles =[
        'Heat',
        'Ion A',
        'Ion B',
        'Ion C',
        'Ion D',
    ]

    quality_cut = df['quality_cut']
    
    source = df['source'].unique()[0]
    num = '{1} ({0}): $\chi^2$ Cut'.format(source, stream)
    
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11.69, 8.27),
                             num=num)

    for suffix, tupl, title in zip(channel_suffix, ax_tuples, ax_titles):
        
        xdata = abs( df['energy_adu_{}'.format(suffix)] )
        ydata = df['chi2_{}'.format(suffix)]
        
        nsamples = xdata.shape[0]
        
        ax = axes[tupl]
        
        ax.plot(xdata, ydata,
                label='All events: {}'.format(nsamples),
                c='red', marker=',', ls='none')
        
        xdata_cut = xdata[quality_cut]
        ydata_cut = ydata[quality_cut]
        
        if nsamples < 1000:
            marker = '.'
        else:
            marker = ','

        ax.plot(xdata_cut, ydata_cut,
                label='Quality events: {}'.format(xdata_cut.shape[0]),
                c='slateblue', marker=marker, ls='none')
    
        ax.legend()
        ax.set_title(title.replace('_', ' '))
        ax.set_xlabel('Energy [ADU]')
        ax.set_ylabel('$\chi^2$')
        ax.set_yscale('log')
        ax.set_xscale('log')
        
        ax.set_xlim(xdata_cut.min()*0.5, ax.get_xlim()[1])
        ax.set_ylim(ydata_cut.min()*0.5, ax.get_ylim()[1])

    fig.text(0.5, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))
        
    fig.delaxes(axes[0,0])    
    fig.tight_layout(rect=(0, 0, 1, 0.98))

    return fig


def histogram_adu(stream, df, bins=1000):
    
        
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
    bulk_cut = df['bulk_cut']
    
    source = df['source'].unique()[0]
    num = '{} ({}) : Quality Cut Histogram'.format(stream, source)

    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(11.69, 8.27),
                             num=num)
    
    for suffix, tupl, label in zip(channel_suffix, ax_tuples, ax_titles):
        
        xdata = df['energy_adu_{}'.format(suffix)]
        
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
        


        xdata_fid = xdata[quality_cut & bulk_cut]
        ax_hist(ax, bin_edges, xdata_fid,
                'Fiducial events', color='limegreen')[0]        
    
### MEANT TO REPRESENT THE 10keV LINE FITTING, not implemented yet
##                if ind in run_tree.chan_signal:
#                if ind == 0: # only for heat channel
#                    
#                    xdata_fit = xdata_qual[(xdata_qual>1000) & (xdata_qual<1400)]
#                    popt = norm.fit(xdata_fit)
#                    
#                    if ana.calibration_peak.cut_type == 'fiducial':
#                        xdata_fit = xdata_fid
#                    
##                    popt = getattr(ana.model.popt, label)
#                    xrange = np.linspace(xdata_fit.min(), xdata_fit.max(), 1000)
##                    pdf = ana.model.dist.pdf(xrange, *popt)
##                    cdf = ana.model.dist.cdf(xrange, *popt)
#                    pdf = norm.pdf(xrange, *popt)
#                    cdf = norm.cdf(xrange, *popt)
#                    normalization = getattr(trig,
#                                            'nsamples_{}'.format(
#                                                    ana.calibration_peak.cut_type
#                                            ))
#                    pdf_norm = pdf * normalization * (bin_edges[1] - bin_edges[0])
#                    
#                    ax.autoscale(False)
#                    ax.plot(xrange, pdf_norm,
#                            ls='--', color='yellow',
#                            label='fit')
#                    
#                    a0.plot(xrange, cdf,
#                            ls='-.', color='yellow',
#                            label='fit')
        
        ax.legend(loc=2)
        ax.set_title(label.replace('_', ' '))
        
    
    fig.text(0.5, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))

### MESSAGE BOX ON THE LEFT, to tell about the calibration, sensitivity etc..
    # msg_list = list()
    # for ind in run_tree.chan_signal:
    #     lab = run_tree.chan_label[ind]
    #     polar, pos10kev, sens_adu, gain_chan, sens_nv, sigma_adu, sigma_ev = optimization_info(ana, ind)
    #     msg_title = r'\underline{'+lab.replace('_', ' ')+'}'
    #     msg_core = (
    #     r'\\ '
    #     r'Position 10keV = {:.0f} ADU \\ '
    #     r'Gain = {} nV/ADU \\ '
    #     r'Sensitivity = {:.2e} ADU/eV = {:.1f} nV/keV \\ '
    #     r'Resolution = {:.2f} ADU = {:.1f} eV '
    #     ).format(pos10kev, gain_chan, sens_adu, sens_nv, sigma_adu, sigma_ev)
        
    #     msg_chan = msg_title+msg_core
    #     msg_list.append(msg_chan)

    # msg = r'\\'.join(msg_list)
    
    # fig.text(0.1, 0.82, msg,
    #          horizontalalignment='left',
    #          verticalalignment='center',
    #          bbox=dict(facecolor='white', alpha=0.5))
    
        
    fig.delaxes(axes[0,0])    
    fig.tight_layout()
        
    return fig


def histogram_ev(stream, df, bins=1000):


    channel_suffix = [
        'heat',
        'ionB',
        'ionD',
        'ion_bulk'
    ]    
    
    ax_titles =[
        'Heat',
        'Ion B',
        'Ion D',
        'Bulk',
    ]
    
    ax_tuples = ((0, 0), (0, 1), (1, 0), (1, 1))

    quality_cut = df['quality_cut']
    bulk_cut = df['bulk_cut']

    source = df['source'].unique()[0]
    num = '{} ({}) : Quality Cut Histogram EV'.format(stream, source)


    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(11.69, 8.27),
                             num=num)
    
    for suffix, tupl, label in zip(channel_suffix, ax_tuples, ax_titles):
        
        xdata = df['energy_{}'.format(suffix)]
        ax = axes[tupl]
        xdata_qual = xdata[quality_cut]
        
        # if etype is trig:
        #     bin_edges = custom_bin_edges(xdata_qual, 
        #                                  getattr(noise.sigma0_ev, label))
    
        bin_edges = np.histogram_bin_edges(xdata[quality_cut], bins=1000)
    
        ax_hist(ax, bin_edges, xdata,
                'All events', color='coral')
        ax_hist(ax, bin_edges, xdata_qual,
                'Quality events', color='slateblue')
        
        xdata_fid = xdata[quality_cut & bulk_cut]
        ax_hist(ax, bin_edges, xdata_fid,
                'Fiducial events', color='limegreen')     
            
        ax.set_xlabel('Enregy [eV]')
        ax.legend(loc=2)
        ax.set_title(label.replace('_', ' '))

    fig.text(0.5, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))
  
    fig.tight_layout()

    return fig


def ion_vs_ion(stream, df):
    
    quality_cut = df['quality_cut']
    bulk_cut = df['bulk_cut']
    
    # initializing pseudo-corner plot
    ax_tuples = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2)]
    ax_discard = [(0, 1), (1, 2), (0, 2)]
    
    # chan_x = np.insert(run_tree.chan_veto, 0, run_tree.chan_collect[1])
    # chan_y = np.append(run_tree.chan_veto, run_tree.chan_collect[0])    
    chan_x = ['ionD', 'ionA', 'ionC']
    chan_y = ['ionA', 'ionC', 'ionB']

    source = df['source'].unique()[0]    
    num = '{} ({}): Ion vs Ion'.format(stream, source)
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(8.27, 8.27),
                             num=num, sharex='col', sharey='row')
    
    # actually plotting the data
    for atupl in ax_tuples:
        
        ax = axes[atupl]
        xind = chan_x[atupl[1]]
        yind = chan_y[atupl[0]]
    
        # energy_x = energy[:, xind]
        # energy_y = energy[:, yind]
        energy_x = df['energy_{}'.format(xind)]
        energy_y = df['energy_{}'.format(yind)]
        
        ax.plot(
                energy_x[quality_cut & bulk_cut], energy_y[quality_cut & bulk_cut],
                ls='none', marker='2', zorder=11, color='limegreen',
                label='Fiducial Events'
        )

        ax.plot(
                energy_x[quality_cut], energy_y[quality_cut],
                ls='none', marker='1', zorder=10, color='slateblue',
                label='Quality Events'
        )
            
        ax.plot(
                energy_x, energy_y,
                ls='none', marker=',', zorder=9, color='coral',
                label='All events'
        )
    
        ### FOR THE DRAWING OF THE BULK CUT
        # if n_sigma != 0:
        #     if xind in run_tree.chan_veto:
        #         lab = run_tree.chan_label[xind]
        #         xamp = n_sigma*getattr(noise.sigma0, lab)
        #         ymin, ymax = energy_y.min(), energy_y.max()
        #         ax.fill_betweenx([ymin, ymax], -xamp, +xamp, color='lavender')
        
        #     if yind in run_tree.chan_veto:
        #         lab = run_tree.chan_label[yind]
        #         yamp = n_sigma*getattr(noise.sigma0, lab)
        #         xmin, xmax = energy_x.min(), energy_x.max()
        #         ax.fill_between([xmin, xmax], -yamp, +yamp, color='lavender',
        #                          label='Fiducial selection ({}$\sigma$)'.format(n_sigma)
        #                         )

            
        custom_autoscale(ax, energy_x[quality_cut], energy_y[quality_cut])
        
        ax.grid(alpha=0.3)
        
        if atupl == (0,0):
            ax.legend(loc='lower left', framealpha=1,
                      bbox_to_anchor=(1.05, 0.05), borderaxespad=0.,
            )
        
        if atupl[0] == 2:
            ax.set_xlabel(
                    'Energy {} [ADU]'.format(xind)
            )
                
        if atupl[1] == 0:
            ax.set_ylabel(
                    'Energy {} [ADU]'.format(yind)
            )
    
    fig.text(0.65, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))
    
    for tupl in ax_discard:
        fig.delaxes(axes[tupl])
    fig.tight_layout()
    fig.subplots_adjust(hspace=.0, wspace=.0)
    
    return fig


def virtual_vs_virtual_ev(stream, df):
    
    quality_cut = df['quality_cut']
    bulk_cut = df['bulk_cut']
    
    
    # initializing pseudo-corner plot
    ax_tuples = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2)]
    ax_discard = [(0, 1), (1, 2), (0, 2)]
    
    
    chan_x = ['heat', 'ionB', 'ionD']
    chan_y = ['ionB', 'ionD', 'ion_bulk']

    source = df['source'].unique()[0]    
    num = '{} ({}): VIRTUAL vs VIRTUAL EV'.format(stream, source)
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(8.27, 8.27),
                             num=num, sharex='col', sharey='row')
    
    # actually plotting the data
    for atupl in ax_tuples:
        
        ax = axes[atupl]
        xlab = chan_x[atupl[1]]
        ylab = chan_y[atupl[0]]
    
        energy_x = df['energy_{}'.format(xlab)]
        energy_y = df['energy_{}'.format(ylab)]
    
        ax.plot(
                energy_x[quality_cut & bulk_cut], energy_y[quality_cut & bulk_cut],
                ls='none', marker='2', zorder=11, color='limegreen',
                label='Fiducial Events'
        )
        
        ax.plot(
                energy_x[quality_cut], energy_y[quality_cut],
                ls='none', marker='1', zorder=10, color='slateblue',
                label='Quality Events'
        )
            
        ax.plot(
                energy_x, energy_y,
                ls='none', marker=',', zorder=9, color='coral',
                label='All events'
        )
            
        custom_autoscale(ax, energy_x[quality_cut], energy_y[quality_cut])
        
        ax.grid(alpha=0.3)
        
        if atupl == (0,0):
            ax.legend(loc='lower left', framealpha=1,
                      bbox_to_anchor=(1.05, 0.05), borderaxespad=0.,
            )
        
        if atupl[0] == 2:
            ax.set_xlabel(
                    'Energy {} [eV]'.format(
                            xlab.replace('_', ' ')
                    )
            )
                
        if atupl[1] == 0:
            ax.set_ylabel(
                    'Energy {} [eV]'.format(
                            ylab.replace('_', ' ')
                    )
            )
    
    fig.text(0.65, 0.98, num,
             horizontalalignment='center',
             verticalalignment='center',
             bbox=dict(facecolor='lime', alpha=0.5))
    
    for tupl in ax_discard:
        fig.delaxes(axes[tupl])
    fig.tight_layout()
    fig.subplots_adjust(hspace=.0, wspace=.0)
    
    return fig

if __name__ == '__main__':
    
    plt.close('all')
    plt.rcParams['text.usetex']=True

    analysis_dir = '/home/misiak/Analysis/neutron_background'
    analysis_data_path = '/'.join([analysis_dir, 'data_analysis.h5'])

    # stream_list = pd.read_hdf(
    #     fine_data_path,
    #     key='df',
    #     columns=['stream',]
    # )['stream'].unique()    
    
    # from tqdm import tqdm
    # for stream in stream_list:
    if True:
        stream = 'tg28l000'
    
        df_analysis = pd.read_hdf(
            analysis_data_path,
            key='df',
            where='stream = "{}"'.format(stream)
        )
        
        # fig_temp = temporal_plot(stream, df_analysis)
        # fig_chi2 = plot_chi2_vs_energy(stream, df_analysis)
        # fig_hist_trig = histogram_adu(stream, df_analysis)
        # fig_hist_trig_ev = histogram_ev(stream, df_analysis)
        # fig_ion = ion_vs_ion(stream, df_analysis)
        fig_virtual = virtual_vs_virtual_ev(stream, df_analysis)
        