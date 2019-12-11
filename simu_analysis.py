#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:27:11 2019

@author: misiak
"""

import numpy as np
import matplotlib.pyplot as plt


#MODE = 'Calibration'
MODE = 'Background'

SAVE_DIR = '/home/misiak/Analysis/fine_root_analysis/fond_neutron/{}/simu'.format(MODE)


npy_dict = {
        'all':'energy_calibrated.npy',
        'quality':'energy_quality.npy',
        'trigger':'energy_trigger.npy',
        'input_all':'input_energy_all.npy',
        'input_quality':'input_energy_quality.npy',
        'input_trigger':'input_energy_trigger.npy',
}

simu_list = {
        'ER':'Flat_Analytical_SimuOnly_0.0000_50.0000_ER',
        'NR':'Flat_Analytical_SimuOnly_0.0000_50.0000_NR',
        '1kev':'Line_Analytical_SimuOnly_1.3000_ER',
        '10kev':'Line_Analytical_SimuOnly_10.3700_ER',
}


D = dict()
for k, v in npy_dict.items():
    
    d= dict()
    for k1, v1 in simu_list.items():
        save_path = '/'.join([SAVE_DIR, v1, k])
        d[k1] = np.load(save_path + '.npy')
    
    D[k] = d


# =============================================================================
# FID CUT
# =============================================================================

std_nb_fid = 1

std_10kev_fid = 1.5

alpha_fid = (std_10kev_fid**2 - std_nb_fid**2)**0.5 / 10.37

def std_fid(collect):
    return ( std_nb_fid**2 + (alpha_fid*collect)**2 )**0.5

def fid_cut(energy, debug=False):
    
    heat, _, A, B, C, D = energy.T
    fid_thresh_a = std_fid(B)
    fid_thresh_c = std_fid(B)
    # HACK THE POSITIVE HEAT IN and the limit to 50keV
    fid_cut = (abs(A)<fid_thresh_a) & (abs(C)<fid_thresh_c) & (heat > 0.025) & (heat < 50)
    
    if debug:
        plt.plot(A, B, ls='none', marker='.', markersize=2, alpha=0.3)
        
        B_graph = np.linspace(B.min(), B.max(), int(1e4))
        
        plt.fill_betweenx(B_graph, -std_fid(B_graph), std_fid(B_graph), alpha=0.3, color='coral')
        
        plt.xlabel('ion A')
        plt.ylabel('ion B')
        plt.grid()

        plt.figure()
        plt.plot(C, D, ls='none', marker='.', markersize=2, alpha=0.3)
        
        D_graph = np.linspace(D.min(), D.max(), int(1e4))
        
        plt.fill_betweenx(D_graph, -std_fid(D_graph), std_fid(D_graph), alpha=0.3, color='coral')
        
        plt.xlabel('ion C')
        plt.ylabel('ion D')
        plt.grid()
        
    return fid_cut


# =============================================================================
# # BAND CUT
# =============================================================================
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

std_nb = 0.254
#std_1kev = 0.272
std_10kev = 0.318

alpha = (std_10kev**2 - std_nb**2)**0.5 / 10.37
def std_collect(ec):
    return ( std_nb**2 + (alpha*ec)**2 )**0.5


er_array = np.linspace(0, 200, int(1e4))
dv=2
ec_array = er_array * (1 + lindhard(er_array)*dv/3) / (1 + dv/3)
ei_array = er_array * lindhard(er_array)

def band_cut(energy):

    heat_fid, heatB,  A_fid, B_fid, C_fid, D_fid = energy.T
    
    collect_fid = (B_fid + D_fid)/2
#    collect = (B + D)/2
    
    
    # gamma cut
    #ec_graph = np.linspace(0, 100, int(1e4))
    #std_graph = std_collect(ec_graph)
    #ax.plot(ec_graph, ec_graph, 
    #        color='deepskyblue',  path_effects=cartoon,
    #        zorder=20)
    #ax.fill_between(ec_graph, ec_graph + 3*std_graph, ec_graph - 3*std_graph,
    #                color='deepskyblue', alpha=0.2,
    #                zorder=-10)
    
    std_fid = std_collect(heat_fid)
    gamma_cut = (collect_fid < (heat_fid + 3*std_fid)) & (collect_fid > (heat_fid - 3*std_fid))
    #ax.plot(heat_fid[gamma_cut], collect_fid[gamma_cut],
    #        label='gamma band',
    #        ls='none', marker='.', alpha=0.7, markersize=10, color='steelblue')
    
    # neutron cut
    #ei_graph = np.interp(ec_graph, ec_array, ei_array)
    #ax.plot(ec_graph, ei_graph, 
    #        color='coral',  path_effects=cartoon,
    #        zorder=20)
    #ax.fill_between(ec_graph,  ei_graph + 3*std_graph, ei_graph - 3*std_graph,
    #                color='coral', alpha=0.2,
    #                zorder=-10)
    
    collect_lindhard = np.interp(heat_fid, ec_array, ei_array)
    neutron_cut = (collect_fid < (collect_lindhard + 3*std_fid)) & (collect_fid > (collect_lindhard - 3*std_fid))
    #ax.plot(heat_fid[neutron_cut], collect_fid[neutron_cut],
    #        label='neutron band',
    #        ls='none', marker='.', alpha=0.7, markersize=10, color='coral')
    
    #ax.legend()
    
    # HO cut
    #ax.plot(ec_graph, np.zeros(ec_graph.shape), 
    #        color='k',  path_effects=cartoon,
    #        zorder=20)
    #ax.fill_between(ec_graph,  3*std_graph, - 3*std_graph,
    #                color='k', alpha=0.2,
    #                zorder=-10)
    
    HO_cut = (collect_fid < (0 + 3*std_fid)) & (collect_fid > (0 - 3*std_fid))
    #ax.plot(heat_fid[HO_cut], collect_fid[HO_cut],
    #        label='HO band',
    #        ls='none', marker='.', alpha=0.7, markersize=10, color='k', zorder=-10)
    #
    #ax.legend()
    
    return gamma_cut, neutron_cut, HO_cut





# =============================================================================
# PLOT
plt.close('all')
# =============================================================================

sim = 'NR'

fid_cut_all = fid_cut(D['all'][sim])
g_cut_all, n_cut_all, ho_cut_all = band_cut(D['all'][sim][fid_cut_all])

fid_cut_trigger= fid_cut(D['trigger'][sim])
g_cut_trigger, n_cut_trigger, ho_cut_trigger = band_cut(D['trigger'][sim][fid_cut_trigger])

#fid_cut = fid_cut(D['quality'][sim], debug=True)
fid_cut = fid_cut(D['quality'][sim], debug=False)
g_cut, n_cut, ho_cut = band_cut(D['quality'][sim][fid_cut])


###
plt.figure('Hist E input')

n_all, bins, p = plt.hist(D['input_all'][sim], bins=100, label='all')

n_trig,_, __ = plt.hist(D['input_trigger'][sim], bins=bins, label='trigger')

n_qual,_, __ = plt.hist(D['input_quality'][sim], bins=bins, label='quality')

n_fid,_, __ = plt.hist(D['input_quality'][sim][fid_cut], bins=bins, label='fid')

n_NR, _, __ = plt.hist(D['input_quality'][sim][fid_cut][n_cut], bins=bins, label='NR')

plt.legend()
plt.xlabel('Input Energy [keV]')


###
plt.figure('Efficacity')

bin_array = bins[:-1] + (bins[1] - bins[0]) / 2

eff_trig = n_trig / n_all
eff_qual = n_qual / n_all
eff_fid = n_fid / n_all
eff_NR = n_NR / n_all

plt.plot(bin_array, eff_trig, label='Trigger + Livetime', ls='steps-mid')

plt.plot(bin_array, eff_qual, label='+ Quality', ls='steps-mid')

plt.plot(bin_array, eff_fid, label='+ FID', ls='steps-mid')

plt.plot(bin_array, eff_NR, label='+ NR band', ls='steps-mid')

plt.grid()
plt.legend()
plt.xlabel('Input Energy [keV]')

