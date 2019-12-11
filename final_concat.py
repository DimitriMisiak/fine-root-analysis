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
        'quality':'energy_quality.npy'
}

mode_list = ['Background', 'Calibration']


with open('/home/misiak/Data/data_run57_neutron/stream_config.json', 'r') as cfg:
    config = json.load(cfg)
            
            
#stream = 'tg27l000'
mode = 'Background'
save_flag = True

for mode in mode_list:

    data_dict = {k:list() for k in npy_dict.keys()}
    for stream in config[mode]:
        save_dir = '/home/misiak/Analysis/fine_root_analysis/fond_neutron/{}'.format(stream)
        
        for k,v in npy_dict.items():
            raw_load = np.load(save_dir+'/'+ v)
    #        data_quality = np.delete(raw_load, 1, axis=1)
            data_dict[k].append(raw_load)
    
    concat_dict = {k:np.concatenate(v) for k,v in data_dict.items()}
    
    for k,v in concat_dict.items():
        save_dir = '/home/misiak/Analysis/fine_root_analysis/fond_neutron/{}'.format(mode)
        os.makedirs(save_dir, exist_ok=True)
        np.save(save_dir + '/' + k, v)




















