#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
@author: misiak

"""

import numpy as np
import matplotlib.pyplot as plt


from plot_addon import LegendTitle, custom_autoscale, ax_hist
from stats_addon import cdf_calc, custom_bin_edges
from scipy.stats import norm

from spec_classes import Analysis_red

def optimization_info(ana, ind):
    
    run_tree = ana.all.run_tree
    
    polar = run_tree.Amp_Chal[0][0] / 0.22e9
    
    lab = run_tree.chan_label[ind]
    pos10kev = getattr(ana.calibration_peak, lab)
    sens_adu = getattr(ana.sensitivity, lab)
    gain_chan = run_tree.Chan_Gain[0][ind]
    sens_nv = sens_adu * gain_chan * 1e3 #nV/keV
    sigma_adu = getattr(ana.all.noise.sigma0, lab)
    sigma_ev = getattr(ana.all.noise.sigma0_ev, lab)

    return polar, pos10kev, sens_adu, gain_chan, sens_nv, sigma_adu, sigma_ev


def temporal_plot(ana):
    
    # general
    run_info = ' '.join([ana.run, ana.detector])
    trig = ana.all.trig
    run_tree = ana.all.run_tree
    
    # cut on data
    cut = np.ones(trig.time.shape, dtype=bool) #no cut
    cut = trig.cut.quality
    
    # data
    energy = trig.filt_decor.Energy_OF
    chi2 = trig.filt_decor.chi2_OF
    offset = trig.raw.Off
    slope = trig.raw.Slope_Ion
    time = trig.time
    
    # Init figure
    num = '{} : Temporal Run Check'.format(run_info)
    fig, axes = plt.subplots(nrows=6, ncols=1, figsize=(11.69, 8.27),
                             sharex=True, num=num)
    
    # heat trig vs time
    ax = axes[0]
    ax.set_ylabel('Energy Heat [ADU]')
    ax.set_yscale('symlog')
    
    ax.plot(
            time[cut], energy[cut, 0],
            label='heat a', zorder=10,
            ls='none', marker='2', mew=0.8,
    )
    ax.autoscale(False)
    ax.plot(
            time, energy[:, 0],
            label='All events',
            ls='none', marker=',', color='silver',
    )
    
    # ion trig vs time
    ax = axes[1]
    ax.set_ylabel('Energy Ion [ADU]')
    ax.set_yscale('symlog')
    
    for i, ind in enumerate(run_tree.chan_ion):
        label = (run_tree.chan_label[ind]).replace('_', ' ')
        ax.plot(
                time[cut], energy[cut, ind],
                label=label, zorder=10,
                ls='none', marker=str(i+1), mew=0.8
        )
    ax.autoscale(False)
    ax.plot(
            time, energy[:, run_tree.chan_ion],
            label='All events',
            ls='none', marker=',', color='silver',
    )
    
    # heat offset vs time
    ax = axes[2]
    ax.plot(
            time[cut], offset[cut, 0],
            label='heat a', zorder=10,
            ls='none', marker='2'
    )
    ax.autoscale(False)
    ax.plot(
            time, offset[:, 0],
            label='All events',
            ls='none', marker=',', c='silver'
    )
    ax.set_ylabel('Offset Heat [ADU]')
    
    # ion offset vs time
    ax = axes[3]
    ax.set_ylabel('Offset Ion [ADU]')
    for i, ind in enumerate(run_tree.chan_ion):
        label = (run_tree.chan_label[ind]).replace('_', ' ')
        ax.plot(
                time[cut], offset[cut, ind],
                label=label, zorder=10,
                ls='none', marker=str(i+1), mew=0.8
        )
    ax.autoscale(False)
    ax.plot(
            time, offset[:, ind],
            label='All events',
            ls='none', marker=',', c='silver'
    )
    
    # ion slope vs time
    ax = axes[4]
    ax.set_ylabel('Slope Ion [ADU/s]')
    for i, ind in enumerate(run_tree.chan_ion):
        label = (run_tree.chan_label[ind]).replace('_', ' ')
        ax.plot(time[cut], slope[cut, ind-2],
                label=label, zorder=10,
                ls='none', marker=str(i+1),
        )
    ax.autoscale(False)
    ax.plot(time, slope[:, ind-2],
            label='All events',
            ls='none', marker=',', c='silver'
    )
    
    # chi2 vs time
    ax = axes[5]
    ax.set_ylabel('$\chi^2$')
    label = 'chi2 heat A'
    ax.plot(time[cut], chi2[cut, 0],
            label=label, zorder=10,
            ls='none', marker=str(i+1),
    )
    ax.autoscale(False)
    ax.plot(time, chi2[:, 0],
            label='All events',
            ls='none', marker=',', c='silver'
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


def plot_chi2_vs_energy(ana, trig_only=False):

    # general
    run_info = ' '.join([ana.run, ana.detector])
    trig = ana.all.trig
    
    run_tree = ana.all.run_tree    
    etypes = [trig,]
    
    if not trig_only:
        noise = ana.all.noise
        etypes = [trig, noise]
    
    fig_list = list()
    for etype in etypes:
    
        if not trig_only:
            if etype is noise:
                energy = etype.filt_decor.Energy_OF_t0
                chi2 = etype.filt_decor.chi2_OF_t0        
            
        if etype is trig:
            energy = etype.filt_decor.Energy_OF
            chi2 = etype.filt_decor.chi2_OF
            
        # chi2 vs Energy plot
        ax_tuples = ((1, 0), (0, 1), (0, 2), (1, 1), (1, 2))       
        data_ind = run_tree.chan_valid
        ax_titles = run_tree.chan_label[data_ind]
        x_datas = (np.abs(energy[:, i]) for i in data_ind)    
        y_datas = (chi2[:, i] for i in data_ind)
    
        num = '{} / {} : Quality Cut Plot'.format(run_info,  etype.name)
        
        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(11.69, 8.27),
                                 num=num)
    
        for tupl, xdata, ydata, title in zip(ax_tuples, x_datas, y_datas, ax_titles):
            
            ax = axes[tupl]
            
            ax.plot(xdata, ydata,
                    label='All events: {}'.format(etype.nsamples),
                    c='red', marker=',', ls='none')
            
            xdata_cut = xdata[etype.cut.quality]
            ydata_cut = ydata[etype.cut.quality]
            
            if etype.nsamples < 1000:
                marker = '.'
            else:
                marker = ','
    
            ax.plot(xdata_cut, ydata_cut,
                    label='Quality events: {}'.format(etype.nsamples_quality),
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

        fig_list.append(fig)

    return fig_list


def histogram_adu(ana, trig_only=False):
    
    # general
    run_info = ' '.join([ana.run, ana.detector])
    trig = ana.all.trig
    run_tree = ana.all.run_tree    
    etypes = [trig,]
    
    if not trig_only:
        noise = ana.all.noise
        etypes = [trig, noise]
        
    fig_list = list()
    for etype in etypes:
    
        if not trig_only:
            if etype is noise:
                energy = etype.filt_decor.Energy_OF_t0
        
        if etype is trig:
            energy = etype.filt_decor.Energy_OF
        
        ax_tuples = ((0, 1), (1, 0), (1, 1), (2, 0), (2, 1))
        data_ind = run_tree.chan_valid
    
        num = '{} / {} : Quality Cut Histogram'.format(run_info,  etype.name)
    
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(11.69, 8.27),
                                 num=num)
        
        for tupl, ind in zip(ax_tuples, data_ind):
            
            xdata = energy[:, ind]
            label = run_tree.chan_label[ind]
            
            ax = axes[tupl]
            xdata_qual = xdata[etype.cut.quality]
            
            if not trig_only:
                if etype is noise:
                    bin_edges = np.histogram_bin_edges(xdata[etype.cut.quality])
            
            if etype is trig:
                try:
                    bin_edges = custom_bin_edges(xdata_qual, 
                                                 getattr(noise.sigma0, label))
                except:
                    bin_edges = np.histogram_bin_edges(xdata[etype.cut.quality])
        
            ax_hist(ax, bin_edges, xdata,
                    'All events', color='coral')
            a0 = ax_hist(ax, bin_edges, xdata_qual,
                    'Quality events', color='slateblue')[0]
            
            if etype is trig:
                
                try:
                    if ana.calibration_peak.cut_type == 'fiducial':
                        xdata_fid = xdata[trig.cut.fiducial]
                        a0 = ax_hist(ax, bin_edges, xdata_fid,
                                'Fiducial events', color='limegreen')[0]        
                except:
                    pass
                
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
            
            if not trig_only:
                if etype is noise:
                    ax.set_yscale('linear')
        
        fig.text(0.5, 0.98, num,
                 horizontalalignment='center',
                 verticalalignment='center',
                 bbox=dict(facecolor='lime', alpha=0.5))
    
        if etype is trig:
            
            try:
                msg_list = list()
                for ind in run_tree.chan_signal:
                    lab = run_tree.chan_label[ind]
                    polar, pos10kev, sens_adu, gain_chan, sens_nv, sigma_adu, sigma_ev = optimization_info(ana, ind)
                    msg_title = r'\underline{'+lab.replace('_', ' ')+'}'
                    msg_core = (
                    r'\\ '
                    r'Position 10keV = {:.0f} ADU \\ '
                    r'Gain = {} nV/ADU \\ '
                    r'Sensitivity = {:.2e} ADU/eV = {:.1f} nV/keV \\ '
                    r'Resolution = {:.2f} ADU = {:.1f} eV '
                    ).format(pos10kev, gain_chan, sens_adu, sens_nv, sigma_adu, sigma_ev)
                    
                    msg_chan = msg_title+msg_core
                    msg_list.append(msg_chan)
            
                msg = r'\\'.join(msg_list)
                
                fig.text(0.1, 0.82, msg,
                         horizontalalignment='left',
                         verticalalignment='center',
                         bbox=dict(facecolor='white', alpha=0.5))
        
            except:
                pass
            
        fig.delaxes(axes[0,0])    
        fig.tight_layout()
        
        fig_list.append(fig)
        
    return fig_list


def histogram_ev(ana):
    
    # general
    run_info = ' '.join([ana.run, ana.detector])
    trig = ana.all.trig
    noise = ana.all.noise
    run_tree = ana.all.run_tree    
    etypes = [trig, noise]
    
    fig_list = list()
    for etype in etypes:
    
        energy = etype.energy_ev
        
        ax_tuples = ((0, 0), (0, 1), (1, 0), (1, 1))
        labels = run_tree.chan_label_virtual
    
        num = '{} / {} : Quality Cut Histogram EV'.format(run_info,  etype.name)
    
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(11.69, 8.27),
                                 num=num)
        
        for tupl, label in zip(ax_tuples, labels):
            xdata = getattr(energy, label)
            ax = axes[tupl]
            xdata_qual = xdata[etype.cut.quality]
            
            if etype is noise:
                bin_edges = np.histogram_bin_edges(xdata[etype.cut.quality])
            
            if etype is trig:
                bin_edges = custom_bin_edges(xdata_qual, 
                                             getattr(noise.sigma0_ev, label))
        
            ax_hist(ax, bin_edges, xdata,
                    'All events', color='coral')
            ax_hist(ax, bin_edges, xdata_qual,
                    'Quality events', color='slateblue')
            
            if etype is trig and ana.calibration_peak.cut_type == 'fiducial':
                xdata_fid = xdata[trig.cut.fiducial]
                ax_hist(ax, bin_edges, xdata_fid,
                        'Fiducial events', color='limegreen')     
                
            ax.set_xlabel('Enregy [eV]')
            ax.legend(loc=2)
            ax.set_title(label.replace('_', ' '))
            
            if etype is noise:
                ax.set_yscale('linear')
        
        fig.text(0.5, 0.98, num,
                 horizontalalignment='center',
                 verticalalignment='center',
                 bbox=dict(facecolor='lime', alpha=0.5))
      
        fig.tight_layout()
        
        fig_list.append(fig)
        
    return fig_list


def ion_vs_ion(ana, n_sigma):
    
    # general
    run_info = ' '.join([ana.run, ana.detector])
    trig = ana.all.trig
    
    run_tree = ana.all.run_tree    
    
    if n_sigma != 0:
        noise = ana.all.noise
    # recovering data
    energy = trig.filt_decor.Energy_OF
    cut_qual = trig.cut.quality
    
    try:
        if ana.calibration_peak.cut_type == 'fiducial':
            cut_fid = trig.cut.fiducial
    except:
        pass
    
    # initializing pseudo-corner plot
    ax_tuples = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2)]
    ax_discard = [(0, 1), (1, 2), (0, 2)]
    
    chan_x = np.insert(run_tree.chan_veto, 0, run_tree.chan_collect[1])
    chan_y = np.append(run_tree.chan_veto, run_tree.chan_collect[0])
    
    num = '{} : Ion vs Ion'.format(run_info)
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(8.27, 8.27),
                             num=num, sharex='col', sharey='row')
    
    # actually plotting the data
    for atupl in ax_tuples:
        
        ax = axes[atupl]
        xind = chan_x[atupl[1]]
        yind = chan_y[atupl[0]]
    
        energy_x = energy[:, xind]
        energy_y = energy[:, yind]
    
        try:
            if ana.calibration_peak.cut_type == 'fiducial':
                ax.plot(
                        energy_x[cut_fid], energy_y[cut_fid],
                        ls='none', marker='2', zorder=11, color='limegreen',
                        label='Fiducial Events'
                )
        except:
            pass
        
        ax.plot(
                energy_x[cut_qual], energy_y[cut_qual],
                ls='none', marker='1', zorder=10, color='slateblue',
                label='Quality Events'
        )
            
        ax.plot(
                energy_x, energy_y,
                ls='none', marker=',', zorder=9, color='coral',
                label='All events'
        )
    
        if n_sigma != 0:
            if xind in run_tree.chan_veto:
                lab = run_tree.chan_label[xind]
                xamp = n_sigma*getattr(noise.sigma0, lab)
                ymin, ymax = energy_y.min(), energy_y.max()
                ax.fill_betweenx([ymin, ymax], -xamp, +xamp, color='lavender')
        
            if yind in run_tree.chan_veto:
                lab = run_tree.chan_label[yind]
                yamp = n_sigma*getattr(noise.sigma0, lab)
                xmin, xmax = energy_x.min(), energy_x.max()
                ax.fill_between([xmin, xmax], -yamp, +yamp, color='lavender',
                                 label='Fiducial selection ({}$\sigma$)'.format(n_sigma)
                                )

            
        custom_autoscale(ax, energy_x[cut_qual], energy_y[cut_qual])
        
        ax.grid(alpha=0.3)
        
        if atupl == (0,0):
            ax.legend(loc='lower left', framealpha=1,
                      bbox_to_anchor=(1.05, 0.05), borderaxespad=0.,
            )
        
        if atupl[0] == 2:
            ax.set_xlabel(
                    'Energy {} [ADU]'.format(
                            run_tree.chan_label[xind].replace('_', ' ')
                    )
            )
                
        if atupl[1] == 0:
            ax.set_ylabel(
                    'Energy {} [ADU]'.format(
                            run_tree.chan_label[yind].replace('_', ' ')
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


def virtual_vs_virtual_ev(ana):

    # general
    run_info = ' '.join([ana.run, ana.detector])
    trig = ana.all.trig
    run_tree = ana.all.run_tree    
        
    if not np.all(np.isin(run_tree.chan_collect, run_tree.chan_signal)):
        return None
        
    # recovering data
    energy = trig.energy_ev
    cut_qual = trig.cut.quality
    
    try:
        if ana.calibration_peak.cut_type == 'fiducial':
            cut_fid = trig.cut.fiducial
    except:
        pass
    
    # initializing pseudo-corner plot
    ax_tuples = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2)]
    ax_discard = [(0, 1), (1, 2), (0, 2)]
    
    
    chan_x = [run_tree.chan_label[ind] for ind in run_tree.chan_collect]
    chan_x.insert(0, 'heat_a')
    chan_y = [run_tree.chan_label[ind] for ind in run_tree.chan_collect]
    chan_y.append('collect')
    
    num = '{} : VIRTUAL vs VIRTUAL EV'.format(run_info)
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(8.27, 8.27),
                             num=num, sharex='col', sharey='row')
    
    # actually plotting the data
    for atupl in ax_tuples:
        
        ax = axes[atupl]
        xlab = chan_x[atupl[1]]
        ylab = chan_y[atupl[0]]
    
        energy_x = getattr(energy, xlab)
        energy_y = getattr(energy, ylab)
    
        try:
            if ana.calibration_peak.cut_type == 'fiducial':
                ax.plot(
                        energy_x[cut_fid], energy_y[cut_fid],
                        ls='none', marker='2', zorder=11, color='limegreen',
                        label='Fiducial Events'
                )
        except:
            pass
        
        ax.plot(
                energy_x[cut_qual], energy_y[cut_qual],
                ls='none', marker='1', zorder=10, color='slateblue',
                label='Quality Events'
        )
            
        ax.plot(
                energy_x, energy_y,
                ls='none', marker=',', zorder=9, color='coral',
                label='All events'
        )
            
        custom_autoscale(ax, energy_x[cut_qual], energy_y[cut_qual])
        
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
    # first command
    plt.close('all')
    plt.rcParams['text.usetex']=True

    run = 'tg25l019'
#    run = 'tg22l000'
    
#    ana = Analysis_red(run, detector='RED80')
    
    ana = Analysis_red(run, detector='RED70',
                         chan_valid=(0, 2, 3, 4, 5),
                         chan_signal=(0,)
                         )   

    fig_temp = temporal_plot(ana)
    fig_chi2_trig, fig_chi2_noise = plot_chi2_vs_energy(ana)
    fig_hist_trig, fig_hist_noise = histogram_adu(ana)
    fig_hist_trig_ev, fig_hist_noise_ev = histogram_ev(ana)
    fig_ion = ion_vs_ion(ana)
    fig_virtual = virtual_vs_virtual_ev(ana)
    
    print(optimization_info(ana, 0))