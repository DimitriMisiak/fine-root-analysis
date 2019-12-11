#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 18:56:20 2019

@author: misiak
"""

import os
import re
import numpy as np
from tqdm import tqdm

stream_list = [
        'tg17l007',
        'tg18l005',
        'tg19l010',
        'tg20l000',
        'tg21l000',
        'tg27l000',
        'tg28l000',
]

simu_list = [
        'Flat_Analytical_SimuOnly_0.0000_50.0000_ER',
        'Flat_Analytical_SimuOnly_0.0000_50.0000_NR',
        'Line_Analytical_SimuOnly_1.3000_ER',
        'Line_Analytical_SimuOnly_10.3700_ER',
]

npy_dict = {
        'all':'energy_calibrated.npy',
        'quality':'energy_quality.npy',
        'trigger':'energy_trigger.npy',
        'input_all':'input_energy_all.npy',
        'input_quality':'input_energy_quality.npy',
        'input_trigger':'input_energy_trigger.npy',
}


def custom_concat(stream, simu):


#    stream = stream_list[0]
#    simu = simu_list[0]
    
    stream_dir = '/home/misiak/Analysis/fine_root_analysis/fond_neutron/{}'.format(stream)
    
    simu_dir = '/'.join([stream_dir, simu + '_tot'])
    os.makedirs(simu_dir, exist_ok=True)
    
    subdir_list = [l for l in os.listdir(stream_dir) if re.match(simu+'_[1-9]', l)]
    
    data_dict = {k:list() for k in npy_dict.keys()}
    
    for subdirname in subdir_list:
        
        dirpath = '/'.join([stream_dir, subdirname])
        for k in npy_dict.keys():
            fname = npy_dict[k]
            filepath = '/'.join([dirpath, fname])
            data_array = np.load(filepath)
            data_dict[k].append(data_array)
    
    concat_dict = {k:np.concatenate(v) for k,v in data_dict.items()}
    
    for k, v in concat_dict.items():
        savepath = '/'.join([simu_dir, npy_dict[k]])
        np.save(savepath, v)


for stream in tqdm(stream_list):
    for simu in simu_list:
        custom_concat(stream, simu)
        


