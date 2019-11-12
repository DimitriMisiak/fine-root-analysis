#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 16:41:10 2019

Exploring the root file containing the data.

@author: misiak
"""

import uproot
import numpy as np
import matplotlib.pyplot as plt
import os

from matplotlib import patheffects as pe
style = [pe.Normal(), pe.withStroke(foreground='k', linewidth=3)]

plt.close('all')

detector_name = 'RED80'

def file_path(run_name, num):
    
    file_name = (
            "ProcessedData_{0}_S{1:02d}_{2}_ChanTrig0.root"
    ).format(run_name, int(num), detector_name)
    
    file_path = DATA_DIR + file_name

    return file_path


# =============================================================================
# HEAT_CHI2_CUT
# =============================================================================
def fun_cut(e):
    cut = 300*np.ones(e.shape)# + (e/100)**2
    return cut


def heat_chi2_cut(run_name, num, chi2_thresh=300):
    """ Return the indexes of the events passing the chi2 cut on the heat
    channel for the given partition.    
    Also return the truth array.
    """
    fp = file_path(run_name, num) 
    root = uproot.open(fp)
 
    tree = root["EventTree_trig_Normal_filt_decor"]
    
    energy_of = tree['Energy_OF'].array()
    energy_chal = energy_of[:, 0].T #keeping ony chalA
    
    chi2_of = tree['chi2_OF'].array()
    chi2_chal = chi2_of[:, 0].T #keeping only chalA
    
    thresh_array = fun_cut(energy_chal)
    
    chi2_cut_cond = chi2_chal < thresh_array

    ind_chi2_cut = np.nonzero(chi2_cut_cond)[0]
    
    return ind_chi2_cut, chi2_cut_cond


def heat_chi2_cut_noise(run_name, num, chi2_thresh=300):
    """ Return the indexes of the events passing the chi2 cut on the heat
    channel for the given partition.    
    Also return the truth array.
    """
    fp = file_path(run_name, num) 
    root = uproot.open(fp)
 
    tree = root["EventTree_noise_Normal_filt_decor"]
    
    energy_of = tree['Energy_OF_t0'].array()
    energy_chal = energy_of[:, 0].T #keeping ony chalA
    
    chi2_of = tree['chi2_OF_t0'].array()
    chi2_chal = chi2_of[:, 0].T #keeping only chalA
    
    thresh_array = fun_cut(energy_chal)
    
    chi2_cut_cond = chi2_chal < thresh_array

    ind_chi2_cut = np.nonzero(chi2_cut_cond)[0]
    
    return ind_chi2_cut, chi2_cut_cond

# =============================================================================
# ION_CHI2_CUT
# =============================================================================

def ion_chi2_cut(run_name, num, chi2_thresh=300):
    """ Return the indexes of the events passing the chi2 cut on the heat
    channel for the given partition.  
    Also return the truth array.
    """
    fp = file_path(run_name, num) 
    root = uproot.open(fp)
 
    tree = root["EventTree_trig_Normal_filt_decor"]
    
    chi2_of = tree['chi2_OF'].array()
    chi2_ion = chi2_of[:, 2:].T #keeping only ion channel
    
    chi2_cut_cond = chi2_ion<chi2_thresh

    cond = chi2_cut_cond[0]
    for cond_aux in chi2_cut_cond[1:]:
        cond = np.logical_and(cond, cond_aux)

    ind_chi2_cut = np.nonzero(cond)[0]
    
    return ind_chi2_cut, cond

def ion_chi2_cut_noise(run_name, num, chi2_thresh=300):
    """ Return the indexes of the events passing the chi2 cut on the heat
    channel for the given partition.  
    Also return the truth array.
    """
    fp = file_path(run_name, num) 
    root = uproot.open(fp)
 
    tree = root["EventTree_noise_Normal_filt_decor"]
    
    chi2_of = tree['chi2_OF_t0'].array()
    chi2_ion = chi2_of[:, 2:].T #keeping only ion channel
    
    chi2_cut_cond = chi2_ion<chi2_thresh

    cond = chi2_cut_cond[0]
    for cond_aux in chi2_cut_cond[1:]:
        cond = np.logical_and(cond, cond_aux)

    ind_chi2_cut = np.nonzero(cond)[0]
    
    return ind_chi2_cut, cond

# =============================================================================
# ION_OFFSET_CUT
# =============================================================================

def ion_offset_cut(run_name, num, off_thresh=14000):
    """ Return the indexes of the events passing the chi2 cut on the heat
    channel for the given partition.    
    Also return the truth array.
    """
    fp = file_path(run_name, num) 
    root = uproot.open(fp)
 
    tree = root["EventTree_trig_Normal_raw"]
    
    offset = tree['Off'].array()
    offset_ion = offset[:, 2:].T #keeping only ion channel
    
    offset_cut_cond = np.abs(offset_ion)<off_thresh

    cond = offset_cut_cond[0]
    for cond_aux in offset_cut_cond[1:]:
        cond = np.logical_and(cond, cond_aux)

    ind_chi2_cut = np.nonzero(cond)[0]
    
    return ind_chi2_cut, cond

def ion_offset_cut_noise(run_name, num, off_thresh=14000):
    """ Return the indexes of the events passing the chi2 cut on the heat
    channel for the given partition.    
    Also return the truth array.
    """
    fp = file_path(run_name, num) 
    root = uproot.open(fp)
 
    tree = root["EventTree_noise_Normal_filt_decor"]
    
    offset = tree['Off'].array()
    offset_ion = offset[:, 2:].T #keeping only ion channel
    
    offset_cut_cond = np.abs(offset_ion)<off_thresh

    cond = offset_cut_cond[0]
    for cond_aux in offset_cut_cond[1:]:
        cond = np.logical_and(cond, cond_aux)

    ind_chi2_cut = np.nonzero(cond)[0]
    
    return ind_chi2_cut, cond


# =============================================================================
# FULL_CUT_NOISE
# =============================================================================

def full_cut_noise(run_name, num):
    """ Return the indexes of the events passing the chi2 cut on the heat
    channel for the given partition.    
    """
    _, cond_heat_chi2_cut = heat_chi2_cut_noise(run_name, num)
    _, cond_ion_chi2_cut = ion_chi2_cut_noise(run_name, num)
    _, cond_ion_off_cut = ion_offset_cut_noise(run_name, num)
    
    cond_full = cond_heat_chi2_cut
    for condi in (cond_ion_chi2_cut, cond_ion_off_cut):
        cond_full = np.logical_and(cond_full, condi)
    
    ind_full_cut = np.nonzero(cond_full)[0]
    
    return ind_full_cut, cond_full

# =============================================================================
# FULL_CUT_EVENT
# =============================================================================

def full_cut_event(run_name, num):
    """ Return the indexes of the events passing the chi2 cut on the heat
    channel for the given partition.    
    """
    _, cond_heat_chi2_cut = heat_chi2_cut(run_name, num)
    _, cond_ion_chi2_cut = ion_chi2_cut(run_name, num)
    _, cond_ion_off_cut = ion_offset_cut(run_name, num)
    
    cond_full = cond_heat_chi2_cut
    for condi in (cond_ion_chi2_cut, cond_ion_off_cut):
        cond_full = np.logical_and(cond_full, condi)
    
    ind_full_cut = np.nonzero(cond_full)[0]
    
    return ind_full_cut, cond_full

save_dir = 'output/run_we_21/RED80'
save_dir_analysis = save_dir + '/analysis'
save_dir_full_cut_event = save_dir + '/full_cut_event'
save_dir_full_cut_noise = save_dir + '/full_cut_noise'

for sd in (save_dir_analysis, save_dir_full_cut_event, save_dir_full_cut_noise):
    try:
        os.makedirs(sd)
    except OSError:
        if not os.path.isdir(sd):
            raise

run_list = (
        'tg19l010',                
        'tg20l000', 
        'tg21l000',
)

energy_of = list()
chi2_of = list()

energy_of_cut = list()
chi2_of_cut = list()

# non decor
energy_of_cut_nd = list()
chi2_of_cut_nd = list()

# offset and slope
offset_list = list()
slope_list = list()

# noise time
noise_stp = list()

# pulse time
pulse_stp = list()

# noise
energy_of_noise = list()
chi2_of_noise = list()

part_in_run_sum = 0
run_duration = 0

for run_name in run_list:    
    
    DATA_DIR = (
            '/home/misiak/Data/data_run57/{0}/{1}/'
    ).format(run_name, detector_name)
    
    part_in_run = len([f for f in os.listdir(DATA_DIR) if "ProcessedData_{}_S".format(run_name) in f])
    
    part_in_run_sum += part_in_run
    
    # =============================================================================
    # ANALYSIS_FULL
    # =============================================================================
    save_dir = save_dir_analysis
    
    ion_label = ('A', 'B', 'C', 'D')
    
    #part_in_run = 3
    
    num_array = np.arange(part_in_run)
    
    for num in num_array:
        fp = file_path(run_name, num) 
        root = uproot.open(fp)
    
        run_tree = root['RunTree_Normal']
        
        tree_raw = root["EventTree_trig_Normal_raw"]
        tree_filt = root["EventTree_trig_Normal_filt"]
        tree_filt_decor = root["EventTree_trig_Normal_filt_decor"]
    
        tree_noise = root["EventTree_noise_Normal_filt_decor"]
    
        # chi2 and energy
        energy_array = tree_filt_decor['Energy_OF'].array()
        chi2_array = tree_filt_decor['chi2_OF'].array()
    
        if energy_array.shape[0] == 0:
            continue
    
        # label non decor = nd
        energy_array_nd = tree_filt['Energy_OF'].array()
        chi2_array_nd = tree_filt['chi2_OF'].array()
    
        energy_of.append(energy_array)
        chi2_of.append(chi2_array)
    
        cut_ind,_ = full_cut_event(run_name, num)
        
        energy_of_cut.append(energy_array[cut_ind])
        chi2_of_cut.append(chi2_array[cut_ind])
        
        energy_of_cut_nd.append(energy_array_nd[cut_ind])
        chi2_of_cut_nd.append(chi2_array_nd[cut_ind])
        
        # offset with cuts
        offset = tree_raw['Off'].array()
        offset_list.append(offset[cut_ind])
        
        # slope with cuts
        slope = tree_raw['Slope_Ion'].array()
        slope_list.append(slope[cut_ind])
        
        # noise with cuts
        energy_array_noise = tree_noise['Energy_OF_t0'].array()
        chi2_array_noise = tree_noise['chi2_OF_t0'].array()
    
        cut_ind_noise,_ = full_cut_noise(run_name, num)
        
        energy_of_noise.append(energy_array_noise[cut_ind_noise])
        chi2_of_noise.append(chi2_array_noise[cut_ind_noise])    
        
        # time stamp with cut
        micro_stp = tree_filt_decor['MicroStp'].array()
        num_part = tree_filt_decor['NumPart'].array()
        freq = run_tree['f_max_heat'].array()[0]
        hour_stp = (micro_stp/3600.) / freq + num_part + run_duration + 1./6
        
        pulse_stp.append(hour_stp[cut_ind])

        # time stamp with cut
        micro_stp = tree_noise['MicroStp'].array()
        num_part = tree_noise['NumPart'].array()
        freq = run_tree['f_max_heat'].array()[0]
        hour_stp = (micro_stp/3600.) / freq + num_part + run_duration + 1./6
        
        noise_stp.append(hour_stp[cut_ind_noise])
        
        # gain
        chan_gain = run_tree['Chan_Gain'].array()[0]
        gain_chal = chan_gain[0]
        gain_ion = chan_gain[2:]
        
        # maintenance
    #    maint_cycle = run_tree['MaintenanceCycle'].array()[0] *1.05 / 3600.
        maint_cycle = run_tree['MaintenanceCycle'].array()[0] / 3600.
        maint_duration = run_tree['MaintenanceDuration'].array()[0] / 3600.

    run_duration = hour_stp[-1]
      
# maintenance time
num_maint_cycle = int(part_in_run_sum // maint_cycle + 1)
    
# data with decor
energy_cut = np.concatenate(energy_of_cut)
chi2_cut = np.concatenate(chi2_of_cut)

stp_cut = np.concatenate(pulse_stp)
noise_stp = np.concatenate(noise_stp)

offset_cut = np.concatenate(offset_list)

slope_cut = np.concatenate(slope_list)
slope_ion = slope_cut.T

# data without decor
energy_cut_nd = np.concatenate(energy_of_cut_nd)
chi2_cut_nd = np.concatenate(chi2_of_cut_nd)

# noise
energy_cut_noise = np.concatenate(energy_of_noise)
chi2_cut_noise = np.concatenate(chi2_of_noise)
    
####
offset_chal = offset_cut[:,0].T
offset_ion = offset_cut[:, 2:].T

energy_cut_chal = energy_cut[:,0].T #keeping only chalA
energy_cut_ion = energy_cut[:,2:].T
chi2_cut_chal = chi2_cut[:, 0].T #keeping only chalA
chi2_cut_ion = chi2_cut[:, 2:].T

energy_cut_chal_nd = energy_cut_nd[:,0].T #keeping only chalA
energy_cut_ion_nd = energy_cut_nd[:,2:].T
chi2_cut_chal_nd = chi2_cut_nd[:, 0].T #keeping only chalA
chi2_cut_ion_nd = chi2_cut_nd[:, 2:].T

energy_chal_noise = energy_cut_noise[:, 0].T
energy_ion_noise = energy_cut_noise[:, 2:].T

chi2_chal_noise = chi2_cut_noise[:, 0].T
chi2_ion_noise = chi2_cut_noise[:, 2:].T

# correcting the sign of the electrodes
energy_cut_conv = np.sum(energy_cut_ion_nd, axis=0)

#volt_config = np.array([[-1, -1, +1, +1]])
v_bias = abs(run_tree['Polar_Ion'].array()[0][1]- run_tree['Polar_Ion'].array()[0][-1])
volt_config = -np.sign(run_tree['Polar_Ion'].array())

energy_cut_ion *= volt_config.T
energy_cut_ion_nd *= volt_config.T

energy_ion_noise *= volt_config.T

  
# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================
#def xmax_histo(bin_array, hist_array):
#    ind_bin_max = np.where(hist_array == hist_array.max())[0]
#    x_bin_max = bin_array[ind_bin_max][0]
#    return x_bin_max

def xmax_histo(bin_array, hist_array, xinf=0):
    condi = bin_array > xinf
    hist = np.where(condi, hist_array, 0)
    ind_bin_max = np.where(hist == hist.max())[0]
    x_bin_max = bin_array[ind_bin_max][0]
    return x_bin_max


###PLOT + CUTS
plt.close('all')

# =============================================================================
# NOISE PLOT
# =============================================================================

#### ION Histogramm check

bin_edges = np.histogram_bin_edges(energy_ion_noise, bins=20)

bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

noise_hist_list = list()
for i in range(4):
    noise_hist, _ = np.histogram(energy_ion_noise[i], bins=bin_edges)
    noise_hist_list.append(noise_hist)

quantiles = np.quantile(energy_ion_noise, (0.16, 0.5, 0.84), axis=1).T

energy_ion_noise_sorted = np.sort(energy_ion_noise, axis=1)
ndim = energy_ion_noise_sorted.shape[-1]
cdf_noise = (np.arange(ndim)+1) / float(ndim)

### PLOT
fig, ax = plt.subplots(nrows=4, figsize=(7,20), sharex=True,
                       num='Noise Ion. energy spectrum')
ax = np.ravel(ax)

sigma_ion = list()

for i,a in enumerate(ax):
    
    lab = ion_label[i]
    
    hist_line, = a.plot(bin_array, noise_hist_list[i], ls='steps-mid',
                        color='slateblue', label='hist')

    a.fill_between(bin_array, noise_hist_list[i], color='deepskyblue',
                    alpha=0.3, step='mid')

    a0 = a.twinx()
    a0.set_ylabel('CDF', color='coral')
    a0.tick_params(axis='y', labelcolor='coral')
    
    cdf_line, = a0.plot(energy_ion_noise_sorted[i], cdf_noise, ls='steps',
                        color='coral', label='cdf', path_effects=style)
    
    a.grid(True)
    a.set_ylabel('Counts Ion {}'.format(lab), color='slateblue')
    a.tick_params(axis='y', labelcolor='k')
    a.set_xlabel('Energy [ADU]')

    quants = quantiles[i]
    
    vline = a.axvline(quants[0], color='k')
    a.axvline(quants[1], color='k', ls='--')
    a.axvline(quants[2], color='k')
    
#    sigma = (-quants[0] + quants[2])/2   
    sigma = np.std(energy_ion_noise[i])
    sigma_ion.append(sigma)
    
    
    a.legend((hist_line, cdf_line, vline, vline), ('hist', 'cdf', '[16th,50th,84th]'),
            title='Ion {0}\n$\sigma_{{ion}}$= {1:.3f} ADU'.format(lab, sigma),
            loc=2, fontsize='small')
    
ax[0].set_title('Ion. Noise with Quality Cuts: {} windows'.format(energy_ion_noise.shape[-1]))

fig.tight_layout()
fig.savefig(save_dir+'/1_ion_noise_hist.pdf')

#### ION B+D Histogramm check

energy_bd_noise = energy_ion_noise[1] + energy_ion_noise[3]

bin_edges = np.histogram_bin_edges(energy_bd_noise, bins=30)

bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

noise_hist_bd, _ = np.histogram(energy_bd_noise, bins=bin_edges)


energy_bd_noise_sorted = np.sort(energy_bd_noise, axis=0)
ndim = energy_bd_noise_sorted.shape[0]
cdf_noise = (np.arange(ndim)+1) / float(ndim)

### PLOT
#fig = plt.figure(figsize=(7,5), num='Ion. B+D Noise energy spectrum')
#ax = fig.subplots()

fig, axis = plt.subplots(nrows=2, figsize=(7,7), num='Ion. B+D and Chal Noise energy spectrum')

ax = axis[0]
ax0 = ax.twinx()

ax.set_title('Ion B+D and Chal Noise with Quality Cuts: {} windows'.format(energy_bd_noise.shape[-1]))

hist_line, = ax.plot(bin_array, noise_hist_bd, ls='steps-mid', color='slateblue')

ax.fill_between(bin_array, noise_hist_bd, color='deepskyblue',
                    alpha=0.3, step='mid')

cdf_line, = ax0.plot(energy_bd_noise_sorted, cdf_noise, ls='steps', color='coral',
                     path_effects=style)

ax.grid(True)
ax.set_ylabel('Counts', color='slateblue')
ax.tick_params(axis='y', labelcolor='slateblue')

ax0.set_ylabel('CDF', color='coral')
ax0.tick_params(axis='y', labelcolor='coral')

ax.set_xlabel('Energy [ADU]')

quants = np.quantile(energy_bd_noise, (0.16, 0.5, 0.84))

vline = ax.axvline(quants[0], color='k')
ax.axvline(quants[1], color='k', ls='--')
ax.axvline(quants[2], color='k')

#sigma = (-quants[0] + quants[2])/2    
sigma = np.std(energy_bd_noise)
sigma_bd = sigma

ax.legend((hist_line, cdf_line, vline), ('hist', 'CDF', '[16th,50th,84th]'),
          title='ION B+D\n$\sigma_{{ion}}$= {:.3f} ADU'.format(sigma),
          loc=2, fontsize='small')

#fig.tight_layout()
#fig.savefig(save_dir+'/ion_bd_noise_hist.pdf')

#### CHAL Histogramm check

popt = np.array([  1.58398155, 842.25988089])

def linear(x,a,b):
    return a*x + b

sens_chal = linear(noise_stp, *popt)/10.37

energy_chal_noise /= sens_chal

bin_edges = np.histogram_bin_edges(energy_chal_noise, bins=30)

bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

noise_hist_chal, _ = np.histogram(energy_chal_noise, bins=bin_edges)


energy_chal_noise_sorted = np.sort(energy_chal_noise, axis=0)
ndim = energy_chal_noise_sorted.shape[0]
cdf_noise = (np.arange(ndim)+1) / float(ndim)

### PLOT
#fig = plt.figure(figsize=(7,5), num='Noise Chal energy spectrum')
#ax = fig.subplots()
ax = axis[1]
ax0 = ax.twinx()

#ax.set_title('ChalA Noise with Quality Cuts: {} windows'.format(energy_chal_noise.shape[-1]))

hist_line, = ax.plot(bin_array, noise_hist_chal, ls='steps-mid', color='slateblue')

ax.fill_between(bin_array, noise_hist_chal, color='deepskyblue',
                    alpha=0.3, step='mid')

cdf_line, = ax0.plot(energy_chal_noise_sorted, cdf_noise, ls='steps', color='coral',
                     path_effects=style)

ax.grid(True)
ax.set_ylabel('Counts', color='slateblue')
ax.tick_params(axis='y', labelcolor='slateblue')

ax0.set_ylabel('CDF', color='coral')
ax0.tick_params(axis='y', labelcolor='coral')

ax.set_xlabel('Energy [ADU]')

quants = np.quantile(energy_chal_noise, (0.16, 0.5, 0.84))

vline = ax.axvline(quants[0], color='k')
ax.axvline(quants[1], color='k', ls='--')
ax.axvline(quants[2], color='k')

sigma = (-quants[0] + quants[2])/2    
sigma = np.std(energy_chal_noise)
sigma_chal = sigma

ax.legend((hist_line, cdf_line, vline), ('hist', 'CDF', '[16th,50th,84th]'),
          title='ChalA\n$\sigma_{{ion}}$= {:.3f} ADU'.format(sigma),
          loc=2, fontsize='small')

fig.tight_layout()
fig.savefig(save_dir+'/2_ion_bd_chal_noise_hist.pdf')

# =============================================================================
# VS Time
# =============================================================================
energy_ion_bd = energy_cut_ion[1] + energy_cut_ion[3]

### PLOT
fig, ax = plt.subplots(nrows=2, figsize=(7,8), sharex=True,
                       num='Ionization B+D and ChalA vs Time')
ax = np.ravel(ax)

ax[0].plot(stp_cut, energy_ion_bd, label='ion B+D', ls='none', marker='+')

ax[1].plot(stp_cut, energy_cut_chal, label='chalA', ls='none', marker='+')

#### evolutive gain
STP = stp_cut[energy_cut_chal>820]
ENE = energy_cut_chal[energy_cut_chal>820]

#ax[1].plot(STP, ENE, label='chalA', ls='none', marker='+')

from scipy.optimize import curve_fit

def linear(x,a,b):
    return a*x + b
    
popt,_ = curve_fit(linear, STP, ENE)

ENE_fit = linear(stp_cut, *popt)

ax[1].plot(stp_cut, ENE_fit, label='fit_gain_chalA', ls='-', color='k')

sens_chal = ENE_fit/10.37
energy_cut_chal /= sens_chal

#### evolutive gain end

ax[0].set_xlabel('Time [hours]')
ax[0].set_ylabel('Energy Ion B+D [ADU]')

ax[1].set_xlabel('Time [hours]')
ax[1].set_ylabel('Energy ChalA [ADU]')

ax[0].set_title('With Quality Cuts: {} events'.format(energy_cut_chal.shape[-1]))

for a in ax:
    a.grid(True)

fig.tight_layout()
fig.savefig(save_dir+'/3_ion_chal_vs_time.pdf')

### PLOT
fig, ax = plt.subplots(nrows=2, figsize=(7,8), sharex=True,
                       num='Offset vs Time')
ax = np.ravel(ax)

for i in range(4):
    lab = ion_label[i]
    ax[0].plot(stp_cut, offset_ion[i], label='Ion {}'.format(lab), ls='-', marker='+')

ax[1].plot(stp_cut, offset_chal, label='ChalA', ls='-', marker='+')

ax[0].set_xlabel('Time [hours]')
ax[0].set_ylabel('Offset Ion [ADU]')

ax[1].set_xlabel('Time [hours]')
ax[1].set_ylabel('Offset Chal [ADU]')

ax[0].set_title(
        ('{0:.2f} hours stream with Quality Cuts: {1} events'
        ).format(run_duration, offset_chal.shape[-1])
)

for a in ax:
    a.legend(loc=0)
    a.grid(True)

fig.tight_layout()
fig.savefig(save_dir+'/4_offset_vs_time.pdf')

### PLOT
fig = plt.figure(figsize=(10,5), num='Slope vs Time')
ax = fig.subplots()

color_list = list()
for i in range(4):
    lab = ion_label[i]
    
#    median = np.median(slope_ion[i])
    
    line, = ax.plot(stp_cut, slope_ion[i], ls='none',
                    marker='+')

    color_list.append(line.get_color())

xlim = ax.get_xlim()
ylim = ax.get_ylim()

beg_list = list()
end_list = list()
beg = 0
for i in range(num_maint_cycle):

    end = beg + maint_cycle/2
    filling = ax.fill_betweenx(ylim, beg, end, color='grey', alpha=0.5)
    
    beg_list.append(beg)
    end_list.append(end)
    
    beg += maint_cycle + maint_duration

stp_good = np.array([])
slope_list = list()

stp_list = list()

for beg, end in zip(beg_list[1:], end_list[:-1]):
    
    cond1 = end <= stp_cut
    cond2 = stp_cut <= beg
    cond_and = np.logical_and(cond1, cond2)
    
    try:
        ind_chunk = np.nonzero(cond_and)[0]
        stp_chunk = stp_cut[ind_chunk]
        slope_chunk = slope_cut[ind_chunk]
        
        ax.plot((stp_chunk[0], stp_chunk[-1]), [0, 0], 'k', lw=3)
        
        stp_list.append((stp_chunk[0], stp_chunk[-1]))
    #    stp_good = np.concatenate((stp_good, stp_chunk))
        slope_list.append(slope_chunk)
    except:
        continue

try:
    slope_good = np.concatenate(slope_list)
    slope_median = np.median(slope_good, axis=0)

    for i in range(4):
        color = color_list[i]
        for j, stp_b in enumerate(stp_list):
            beg, end = stp_b
            if j == 0:
                ax.plot([beg, end], [slope_median[i], slope_median[i]], 
                        color=color, path_effects=style, zorder=10,
                        label='$I_{{Ion{0}}}$={1:.2e} ADU/s'.format(lab, slope_median[i]))    
                continue
            ax.plot([beg, end], [slope_median[i], slope_median[i]], 
            color=color, path_effects=style, zorder=10)   
    
    filling.set_label('Maintenance')
    
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    
    ax.set_xlabel('Time [hours]')
    ax.set_ylabel('Slope [ADU/s]')
    
    ax.set_title('With Quality Cuts: {} events'.format(slope_ion.shape[-1]))
    
    ax.legend(title='Median Leakage Currents:', loc=0)
    
    ax.grid(True)
    
    fig.tight_layout()
    fig.savefig(save_dir+'/5_ion_slope_vs_time.pdf')
    
except:
    pass

# =============================================================================
# CRYSTAL EVENT CUT PLOT + FIDUCIAL CUT
# =============================================================================
ion_config = ((1,3), (2,3), (1,0), (2,0))

# fiducial cut
sigma_0 = 1.241
sigma_calib = 3.7
adu_calib = 52.6
alpha = ((sigma_calib**2-sigma_0**2)/adu_calib**2)**0.5

def evol_sigma(sigma, e_adu):
    sigma_evo = (sigma**2 + (alpha*e_adu)**2)**0.5
    return sigma_evo

e_cut_A = 2*evol_sigma(sigma_ion[0], energy_cut_ion[0])
e_cut_C = 2*evol_sigma(sigma_ion[2], energy_cut_ion[2])

def fiducial_cut(e_ion):
    
    cond_a_sup = e_ion[0] < e_cut_A
    cond_a_inf = e_ion[0] > - e_cut_A

    cond_c_sup = e_ion[2] < e_cut_C
    cond_c_inf = e_ion[2] >  - e_cut_C
    
    cond_fid = cond_a_sup
    for condi in (cond_a_inf, cond_c_inf, cond_c_sup):
        cond_fid = np.logical_and(cond_fid, condi)
    
    ind_fid_cut = np.nonzero(cond_fid)[0]
    
    return ind_fid_cut, cond_fid

fid_cut, _ = fiducial_cut(energy_cut_ion)

energy_fid_ion = energy_cut_ion[:, fid_cut]
energy_fid_chal = energy_cut_chal[fid_cut]

### PLOT
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10,10), num='Ion vs Ion',
                       sharex=True, sharey=True)
ax = np.ravel(ax)

for i, config in enumerate(ion_config):
    
    ionx, iony = config

    ax[i].plot(energy_cut_ion_nd[ionx], energy_cut_ion_nd[iony], label='all events (non decor)',
               color='grey', ls='none', marker='.', alpha=0.3)
    
    ax[i].plot(energy_cut_ion[ionx], energy_cut_ion[iony], 
               color='r', ls='none', marker='+',
               label='all events = {}'.format(energy_cut_ion.shape[-1]))

    ax[i].plot(energy_fid_ion[ionx], energy_fid_ion[iony], 
               color='slateblue', ls='none', marker='+',
               label='fiducial events = {}'.format(energy_fid_ion.shape[-1]))

    ax[i].set_ylabel('Ion {} [ADU]'.format(ion_label[iony]))
    ax[i].set_xlabel('Ion {} [ADU]'.format(ion_label[ionx]))
    ax[i].grid(True)

try:
    xylim = (energy_cut_ion.min()*1.2, energy_cut_ion.max()*1.2)
    
    vline = [y*1.1 for y in xylim]
    hline = [x*1.1 for x in xylim]
#    e_inf = [e_cut_inf, e_cut_inf]
#    e_sup = [e_cut_sup, e_cut_sup]
#    
    vline = np.linspace(*xylim, 100)
    e_A = 2*evol_sigma(sigma_ion[0], vline)
    e_C = 2*evol_sigma(sigma_ion[2], vline)
    
    ax[2].fill_between(vline, -e_A, e_A, color='lime', alpha=0.5, 
                       label='fiducial cut on A')
    ax[3].fill_between(vline, -e_A, e_A, color='lime', alpha=0.5)

    ax[1].fill_betweenx(vline, -e_C, e_C, color='lime', alpha=0.5, 
                       label='fid cut on C')
    ax[3].fill_betweenx(vline, -e_C, e_C, color='lime', alpha=0.5)
    
#    ax[1].fill_betweenx(vline, e_inf, e_sup, color='lime', alpha=0.5)
#    ax[3].fill_betweenx(vline, e_inf, e_sup, color='lime', alpha=0.5)
    
    ax[2].legend(title='With quality cuts:')
    
    ax[0].set_title('Ion vs Ion with Quality Cuts')
    ax[0].set_xlim(*xylim)
    ax[0].set_ylim(*xylim)
    
    fig.tight_layout()
    fig.savefig(save_dir+'/6_ion_vs_ion.pdf')

except:
    pass
# =============================================================================
# Ionization Histogramm
# =============================================================================

#### Histogramm check
bin_edges = np.histogram_bin_edges((energy_cut_ion), bins=100)

bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

cut_hist_list = list()
fid_hist_list = list()
for i in range(4):
    cut_hist, _ = np.histogram(energy_cut_ion[i], bins=bin_edges)
    cut_hist_list.append(cut_hist)
    
    fid_hist, _ = np.histogram(energy_fid_ion[i], bins=bin_edges)
    fid_hist_list.append(fid_hist)

energy_cut_ion_sorted = np.sort(energy_cut_ion, axis=1)
ndim = energy_cut_ion_sorted.shape[-1]
cdf_cut = (np.arange(ndim)+1) / float(ndim)

energy_fid_ion_sorted = np.sort(energy_fid_ion, axis=1)
ndim = energy_fid_ion_sorted.shape[-1]
cdf_fid = (np.arange(ndim)+1) / float(ndim)

#### FITTING THE CDF
from scipy.stats import norm
from scipy.optimize import curve_fit

loc_list = list()
scale_list = list()

#plt.figure('testy')
for i in range(4):
    
    ind = np.nonzero(energy_cut_ion[i]>40)[0]
    loc, scale = norm.fit(energy_cut_ion[i][ind])
    loc_list.append(loc)
    scale_list.append(scale)
#
#    plt.hist(energy_cut_ion[i][ind])
#    x_data = np.linspace(0, 100, 1000)
#    mod = norm.cdf(x_data, loc, scale)
#    plt.plot(x_data, mod)


### PLOT
fig, ax = plt.subplots(nrows=4, figsize=(7,10), sharex=True,
                       num='Ion. energy spectrum')
ax = np.ravel(ax)

for i,a in enumerate(ax):
    
    lab = ion_label[i]
    
    hist_line, = a.plot(bin_array, cut_hist_list[i], ls='steps-mid', color='r')

    a.fill_between(bin_array, cut_hist_list[i], color='coral',
                   alpha=0.3, step='mid')

    hist_line_2, = a.plot(bin_array, fid_hist_list[i], ls='steps-mid', color='slateblue')

    a.fill_between(bin_array, fid_hist_list[i], color='deepskyblue',
                   alpha=0.3, step='mid')

    a0 = a.twinx()
    a0.set_ylabel('CDF', color='grey')
    a0.tick_params(axis='y', labelcolor='grey')
    
    cdf_line, = a0.plot(energy_cut_ion_sorted[i], cdf_cut, ls='steps', color='coral',
                        path_effects=style)


    loc = loc_list[i]
    scale = scale_list[i]
    fit_line, = a0.plot(energy_cut_ion_sorted[i],
                        norm.cdf(energy_cut_ion_sorted[i], loc, scale),
                        color='k'
                        )

    cdf_line_2, = a0.plot(energy_fid_ion_sorted[i], cdf_fid, ls='steps', color='deepskyblue',
                        path_effects=style)
    
    a.grid(True)
    a.set_ylabel('Counts Ion {}'.format(lab), color='k')
    a.tick_params(axis='y', labelcolor='k')
    a.set_xlabel('Energy [ADU]')
#    a.set_yscale('log')

    sens = xmax_histo(bin_array, fid_hist_list[i])
    a.axvline(sens, color='k', ls='-.', lw=0.5)

    a0.legend((hist_line, cdf_line, hist_line_2, cdf_line_2, fit_line), 
              ('Hist quality evts', 'Cdf Quality evts', 'Hist FID evts', 'Cdf FID evts', 'Fit:$\mu={:.1f}$, $\sigma={:.1f}$'.format(loc,scale)),
              title=(
                      'Ion {0}\n$S_{{V,ion}}$= {1:.3f} ADU/keV\nGain={2}nV/ADU\nSv={3:.1f}nV/keV'
              ).format(lab, sens/10.37, gain_ion[i], sens/10.37 * gain_ion[i]),
              loc=0, fontsize='small'
    )

ax[0].set_title(
        ('Ion Channels: {} Quality Events, {} Fiducial Events'
         ).format(energy_cut_ion.shape[-1], energy_fid_ion.shape[-1])
)

fig.tight_layout()
fig.savefig(save_dir+'/7_ion_energy_spectrum.pdf')

# =============================================================================
# Heat Histogramm 
# =============================================================================

### Histogramm check

bin_edges = np.histogram_bin_edges(energy_cut_chal, bins=250)

bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

hist_cut_chal, _ = np.histogram(energy_cut_chal, bins=bin_edges)
energy_cut_chal_sorted = np.sort(energy_cut_chal, axis=0)
ndim_cut = energy_cut_chal_sorted.shape[0]
cdf_cut = (np.arange(ndim_cut)+1) / float(ndim_cut)

hist_fid_chal, _ = np.histogram(energy_fid_chal, bins=bin_edges)
energy_fid_chal_sorted = np.sort(energy_fid_chal, axis=0)
ndim_fid = energy_fid_chal_sorted.shape[0]
cdf_fid = (np.arange(ndim_fid)+1) / float(ndim_fid)

### FITTING CDF
ind = np.nonzero(energy_fid_chal>9)[0]
loc_chal, scale_chal = norm.fit(energy_fid_chal[ind])

#### PLOT
#fig = plt.figure(figsize=(7,5), num='Chal energy spectrum')
#ax = fig.subplots()


fig, axis = plt.subplots(nrows=2, figsize=(10,7), num='Ion. B+D and Chal energy spectrum')


### CHAL HISTO
ax = axis[1]
ax0 = ax.twinx()

ax.set_title('Chal A')

hist_line, = ax.plot(bin_array, hist_cut_chal, ls='steps-mid', color='r')

ax.fill_between(bin_array, hist_cut_chal, color='coral', alpha=0.3, step='mid')

hist_line_2, = ax.plot(bin_array, hist_fid_chal, ls='steps-mid', color='slateblue')

ax.fill_between(bin_array, hist_fid_chal, color='deepskyblue', alpha=0.3, step='mid')

cdf_line, = ax0.plot(energy_cut_chal_sorted, cdf_cut, ls='steps', color='coral',
                     path_effects=style)

cdf_line_2, = ax0.plot(energy_fid_chal_sorted, cdf_fid, ls='steps', color='deepskyblue',
                       path_effects=style)

fit_line, = ax0.plot(energy_fid_chal_sorted,
                    norm.cdf(energy_fid_chal_sorted, loc_chal, scale_chal),
                    color='k'
                    )

ax.grid(True)
ax.set_ylabel('Counts', color='k')
ax.tick_params(axis='y', labelcolor='k')

ax0.set_ylabel('CDF', color='grey')
ax0.tick_params(axis='y', labelcolor='grey')

ax.set_xlabel('Energy [ADU]')
#ax.set_xlim(energy_cut_chal_sorted.min(), energy_cut_chal_sorted.max())

sens = xmax_histo(bin_array, hist_fid_chal, xinf=8)
ax.axvline(sens, color='k', ls='-.', lw=0.5)

ax.legend((hist_line, cdf_line, hist_line_2, cdf_line_2, fit_line), 
          ('Hist quality evts', 'Cdf Quality evts', 'Hist FID evts', 'Cdf FID evts', 'Fit:$\mu={:.1f}$, $\sigma={:.1f}$'.format(loc_chal,scale_chal)),
          title='ChalA\n$S_{{V,chal}}$= {:.3f} ADU/keV\nGain={}nV/ADU\nSv={:.1f}nV/keV'.format(sens, gain_chal, sens/10.37 * gain_chal),
          loc=2, fontsize='small')

#fig.tight_layout()
#fig.savefig(save_dir+'/chal_energy_spectrum.pdf')

# =============================================================================
# HISTOGRAMM B+D
# =============================================================================
energy_ion_bd_fid = energy_fid_ion[1] + energy_fid_ion[3]

### Histogramm check

bin_edges = np.histogram_bin_edges(energy_ion_bd , bins=200)

bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

hist_ion_bd, _ = np.histogram(energy_ion_bd, bins=bin_edges)
hist_ion_bd_fid, _ = np.histogram(energy_ion_bd_fid, bins=bin_edges)

energy_ion_bd_sorted = np.sort(energy_ion_bd, axis=0)
ndim_bd = energy_ion_bd_sorted.shape[0]
cdf_bd = (np.arange(ndim_bd)+1) / float(ndim_bd)

energy_ion_bd_sorted_fid = np.sort(energy_ion_bd_fid, axis=0)
ndim_bd_fid = energy_ion_bd_sorted_fid.shape[0]
cdf_bd_fid = (np.arange(ndim_bd_fid)+1) / float(ndim_bd_fid)

### FITTING CDF
### FITTING CDF
ind = np.nonzero(energy_ion_bd_fid>80)[0]
loc_ion, scale_ion = norm.fit(energy_ion_bd_fid[ind])

### PLOT
#fig = plt.figure(figsize=(7,5), num='Ion B+D energy spectrum')
#ax = fig.subplots()

ax = axis[0]
ax0 = ax.twinx()

ax.set_title(
        ('Ion. B+D: {} Quality Events, {} Fiducial Events'
         ).format(energy_cut_chal.shape[-1], energy_fid_chal.shape[-1])
)

hist_line, = ax.plot(bin_array, hist_ion_bd, ls='steps-mid', color='r')
ax.fill_between(bin_array, hist_ion_bd, color='coral', alpha=0.3, step='mid')

hist_line_2, = ax.plot(bin_array, hist_ion_bd_fid, ls='steps-mid', color='slateblue')
ax.fill_between(bin_array, hist_ion_bd_fid, color='deepskyblue', alpha=0.3, step='mid')

cdf_line, = ax0.plot(energy_ion_bd_sorted, cdf_bd, ls='steps', color='coral',
                     path_effects=style)

cdf_line, = ax0.plot(energy_ion_bd_sorted_fid, cdf_bd_fid, ls='steps', color='deepskyblue',
                     path_effects=style)

fit_line, = ax0.plot(energy_ion_bd_sorted_fid,
                    norm.cdf(energy_ion_bd_sorted_fid, loc_ion, scale_ion),
                    color='k'
                    )

ax.grid(True)
ax.set_ylabel('Counts', color='k')
ax.tick_params(axis='y', labelcolor='k')

ax0.set_ylabel('CDF', color='grey')
ax0.tick_params(axis='y', labelcolor='grey')

ax.set_xlabel('Energy [ADU]')

sens_bd = xmax_histo(bin_array, hist_ion_bd_fid)
ax.axvline(sens_bd, color='k', ls='-.', lw=0.5)

ax.legend((hist_line, cdf_line, hist_line_2, cdf_line_2, fit_line), 
          ('Hist quality evts', 'Cdf Quality evts', 'Hist FID evts', 'Cdf FID evts', 'Fit:$\mu={:.1f}$, $\sigma={:.1f}$'.format(loc_ion,scale_ion)),
          title='Ion. B+D\n$S_{{V,ion}}$= {:.3f} ADU/keV\nGain={}nV/ADU\nSv={:.1f}nV/keV'.format(sens_bd/10.37, gain_ion[0], sens_bd/10.37 * gain_ion[0]),
          loc=2, fontsize='small')

fig.tight_layout()
fig.savefig(save_dir+'/8_ion_bd_chal_energy_spectrum.pdf')

# =============================================================================
# ION B+D vs Chal
# =============================================================================
#%%
#plt.close('all')
alpha_chal = ((scale_chal**2-sigma_chal**2)/loc_chal**2)**0.5

def evol_sigma_chal(energy):
    sigma_evo = (sigma_chal**2 + (alpha_chal*energy)**2)**0.5
    return sigma_evo

sens_ion = (loc_ion/10.37)
sigma_bd /= sens_ion
scale_ion /= sens_ion
alpha_bd = ((scale_ion**2-sigma_bd**2)/10.37**2)**0.5

def evol_sigma_bd(energy):
    sigma_evo = (sigma_bd**2 + (alpha_bd*energy)**2)**0.5
    return sigma_evo

from scipy.stats import uniform

er_array = uniform.rvs(scale=25, size=10000)

Q_er = 1
Q_nr = 0.16* er_array**0.18

loc_ei_er = er_array * Q_er
loc_ei_nr = er_array * Q_nr

loc_ec_er = er_array * (1 + Q_er * v_bias / 3)
loc_ec_nr = er_array * (1 + Q_nr * v_bias / 3)

sc_array = evol_sigma_chal(loc_ec_er)
si_array = evol_sigma_bd(loc_ei_er)

ei_er_array = norm.rvs(loc=loc_ei_er, scale=2*si_array)
ei_nr_array = norm.rvs(loc=loc_ei_nr, scale=2*si_array)

ec_er_array = norm.rvs(loc=loc_ec_er, scale=2*sc_array)
ec_nr_array = norm.rvs(loc=loc_ec_nr, scale=2*sc_array)

### PLOT
fig = plt.figure(figsize=(7,5), num='Ionization B+D vs ChalA')
ax = fig.subplots()

ax.plot(energy_cut_chal, energy_ion_bd/sens_ion, label='quality events = {}'.format(energy_cut_chal.shape[-1]),
        color='r', ls='none', marker='+', alpha=0.3, zorder=10)

ax.plot(energy_fid_chal, energy_ion_bd_fid/sens_ion, label='fiducial events = {}'.format(energy_ion_bd_fid.shape[-1]),
        color='slateblue', ls='none', marker='+', zorder=10)

#ax.plot(ec_er_array, ei_er_array, label='ER Monte-Carlo',
#        color='lime', ls='none', marker='o', alpha=0.3)
#
#ax.plot(ec_nr_array, ei_nr_array, label='NR Monte-Carlo',
#        color='crimson', ls='none', marker='o', alpha=0.3)

#ax.plot(ec_er_array, ei_er_array, label='ER Monte-Carlo',
#        color='lime', ls='none', marker='o', alpha=0.3)
#
#ax.plot(ec_nr_array, ei_nr_array, label='NR Monte-Carlo',
#        color='crimson', ls='none', marker='o', alpha=0.3)


_, xsup = ax.get_xlim()
ylim = ax.get_ylim()

ax.set_title('Ionization B+D vs ChalA')

ax.axhline(sens_bd, color='k', ls='-.', lw=0.5)
ax.axvline(sens, color='k', ls='-.', lw=0.5)

#ax.plot([0,sens*2], [0, sens_bd*2], color='coral', alpha=0.3, lw=10, label='ER band (visual guide)')
#ax.plot([0,sens*2], [0, 0.3*2*sens_bd], color='deepskyblue', alpha=0.3, lw=10, label='NR band (Q=0.3, visual guide)')

ax.set_xlim(-xsup/20, xsup)
ax.set_ylim(*ylim)

ax.set_xlabel('Energy ChalA [keV]')
ax.set_ylabel('Energy Ion B+D [keV]')
#ax.legend(title='$S_{{V,ion}}$= {:.3f} ADU/keV\n$S_{{V,chal}}$= {:.3f} ADU/keV\nSv={:.1f}nV/keV'.format(sens_bd/10.37, sens/10.37, sens_bd/10.37 * sens/10.37),
#          fontsize='small', loc=0)
ax.legend(loc=0)

ax.grid(True)

fig.tight_layout()
fig.savefig(save_dir+'/9_ion_vs_chal.pdf')


#%%

## =============================================================================
## Charge Conservation
## =============================================================================
### Histogramm check

energy_fid_conv = energy_cut_conv[fid_cut]

bin_edges = np.histogram_bin_edges(energy_cut_conv, bins=30)

bin_array = bin_edges[1:] - (bin_edges[1]-bin_edges[0])/2

hist_cut_conv, _ = np.histogram(energy_cut_conv, bins=bin_edges)
energy_cut_conv_sorted = np.sort(energy_cut_conv, axis=0)
ndim_cut = energy_cut_conv_sorted.shape[0]
cdf_cut = (np.arange(ndim_cut)+1) / float(ndim_cut)

hist_fid_conv, _ = np.histogram(energy_fid_conv, bins=bin_edges)
energy_fid_conv_sorted = np.sort(energy_fid_conv, axis=0)
ndim_fid = energy_fid_conv_sorted.shape[0]
cdf_fid = (np.arange(ndim_fid)+1) / float(ndim_fid)

#### PLOT
fig = plt.figure(figsize=(7,5), num='Charge conservation: Ion A+B+C+D energy spectrum')
ax = fig.subplots()

### conv HISTO
#ax = axis[1]
ax0 = ax.twinx()

ax.set_title('Charge conservation: Ion A+B+C+D energy spectrum')

hist_line, = ax.plot(bin_array, hist_cut_conv, ls='steps-mid', color='r')

ax.fill_between(bin_array, hist_cut_conv, color='coral', alpha=0.3, step='mid')

hist_line_2, = ax.plot(bin_array, hist_fid_conv, ls='steps-mid', color='slateblue')

ax.fill_between(bin_array, hist_fid_conv, color='deepskyblue', alpha=0.3, step='mid')

cdf_line, = ax0.plot(energy_cut_conv_sorted, cdf_cut, ls='steps', color='coral',
                     path_effects=style)

cdf_line_2, = ax0.plot(energy_fid_conv_sorted, cdf_fid, ls='steps', color='deepskyblue',
                       path_effects=style)

ax.grid(True)
ax.set_ylabel('Counts', color='k')
ax.tick_params(axis='y', labelcolor='k')

ax0.set_ylabel('CDF', color='grey')
ax0.tick_params(axis='y', labelcolor='grey')

ax.set_xlabel('Energy [ADU]')
#ax.set_xlim(energy_cut_conv_sorted.min(), energy_cut_conv_sorted.max())

#sens = xmax_histo(bin_array, hist_fid_conv)
#ax.axvline(sens, color='k', ls='-.', lw=0.5)

cut_mean = np.mean(energy_cut_conv)
cut_sigma = np.std(energy_cut_conv)

fid_mean = np.mean(energy_fid_conv)
fid_sigma = np.std(energy_fid_conv)

ax.legend((hist_line, cdf_line, hist_line_2, cdf_line_2), 
          ('Hist Quality evts', 'Cdf Quality evts', 'Hist FID evts', 'Cdf FID evts'),
          title='Qual Mean = {0:.2f} ADU\nQual Sigma = {1:.2f} ADU\nFID Mean = {2:.2f} ADU\nFID Sigma = {3:.2f} ADU'.format(cut_mean, cut_sigma, fid_mean, fid_sigma),
          loc=0, fontsize='small')

fig.tight_layout()
fig.savefig(save_dir+'/91_charge_conservation.pdf')                                                                                             

# =============================================================================
# Quenching plot
# =============================================================================

# calib coeff ADU/keV
calib_chal = 1 / ( (10.37) * (1 + v_bias / 3) )
calib_ion = sens_bd/10.37
   

ER_cut = energy_cut_chal/calib_chal - (energy_ion_bd/calib_ion) * v_bias / 3
Q_cut = (energy_ion_bd/calib_ion)/ER_cut

ER_fid = energy_fid_chal/calib_chal - (energy_ion_bd_fid/calib_ion) * v_bias / 3
Q_fid = (energy_ion_bd_fid/calib_ion)/ER_fid    

#    ER_cut = energy_cut_chal/calib_chal - (energy_ion_bd/calib_ion) * v_bias * 1.6e-19 / 0.003
#    Q_cut = (energy_ion_bd/calib_ion)/ER_cut
#    
#    ER_fid = energy_fid_chal/calib_chal - (energy_ion_bd_fid/calib_ion) * v_bias * 1.6e-19 / 0.003
#    Q_fid = (energy_ion_bd_fid/calib_ion)/ER_fid

ER_array = np.linspace(0, 15, 1000)
Q_array = 0.16* ER_array**0.18
### PLOT
#    fig = plt.figure(figsize=(7,5), num='Quenching Plot')
#    ax = fig.subplots()

fig, axis = plt.subplots(ncols=2, figsize=(18,7), num='Quenching Plot')

for ax in axis:

    ax.plot(ER_cut, Q_cut, label='quality events = {}'.format(energy_cut_chal.shape[-1]),
            color='r', ls='none', marker='+', alpha=0.3, zorder=10)
    
    ax.plot(ER_fid, Q_fid, label='fid events = {}'.format(energy_fid_chal.shape[-1]),
            color='slateblue', ls='none', marker='+', zorder=10)

    ax.plot(ER_array, Q_array, label='$0.16*E_R^{0.18}$',
            color='lime', alpha=0.3, lw=10)
  
    xinf, xsup = ax.get_xlim()
#    ylim = ax.get_ylim()
    
    ax.set_title('Quenching Plot Q($E_R$) GOOD')
    
    #ax.plot([0,sens*2], [0, sens_bd*2], color='coral', alpha=0.3, lw=10, label='ER band (visual guide)')
    #ax.plot([0,sens*2], [0, 0.3*2*sens_bd], color='deepskyblue', alpha=0.3, lw=10, label='NR band (Q=0.3, visual guide)')
    
#    ax.set_xlim(xinf, 12)
#    ax.set_ylim(-0.5, 1.5)
    
    ax.set_xlabel('Recoil Energy $E_R = E_C - E_I q \Delta V / 3eV$ [keV]')
    ax.set_ylabel('Quenching coeff. $Q=E_I/E_R$')
    ax.legend()
    ax.grid(True)
    
axis[1].set_xscale('log')  
axis[1].set_xlim(0.5, 12)
fig.tight_layout()
fig.savefig(save_dir+'/quenching_plot.pdf')

