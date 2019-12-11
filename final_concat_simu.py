#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 13:36:56 2019

@author: misiak
"""



import numpy as np
import json
import os

npy_dict = {
        'all':'energy_calibrated.npy',
        'quality':'energy_quality.npy',
        'trigger':'energy_trigger.npy',
        'input_all':'input_energy_all.npy',
        'input_quality':'input_energy_quality.npy',
        'input_trigger':'input_energy_trigger.npy',
}

simu_list = [
        'Flat_Analytical_SimuOnly_0.0000_50.0000_ER',
        'Flat_Analytical_SimuOnly_0.0000_50.0000_NR',
        'Line_Analytical_SimuOnly_1.3000_ER',
        'Line_Analytical_SimuOnly_10.3700_ER',
]

mode_list = ['Background', 'Calibration']


with open('/home/misiak/Data/data_run57_neutron/stream_config.json', 'r') as cfg:
    config = json.load(cfg)
            
            
#stream = 'tg27l000'
mode = 'Background'
simu = simu_list[0]
save_flag = True

for mode in mode_list:
    for simu in simu_list:


        
        data_dict = {k:list() for k in npy_dict.keys()}
        for stream in config[mode]:
            save_dir = '/home/misiak/Analysis/fine_root_analysis/fond_neutron/{}/{}_tot'.format(stream, simu)
            
            for k,v in npy_dict.items():
                raw_load = np.load(save_dir+'/'+ v)
        #        data_quality = np.delete(raw_load, 1, axis=1)
                data_dict[k].append(raw_load)
        
        concat_dict = {k:np.concatenate(v) for k,v in data_dict.items()}
        
        for k,v in concat_dict.items():
            save_dir = '/home/misiak/Analysis/fine_root_analysis/fond_neutron/{}/simu/{}'.format(mode, simu)
            os.makedirs(save_dir, exist_ok=True)
            np.save(save_dir + '/' + k, v)




















