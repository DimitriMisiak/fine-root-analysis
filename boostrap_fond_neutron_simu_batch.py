#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 18:56:20 2019

@author: misiak
"""

import os

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

simu_list_1 = list()
for simu in simu_list:
    simu_1 = [simu+'_'+str(x) for x in range(1,10)]
    [simu_list_1.append(x) for x in simu_1]


for stream in stream_list:
    for simu in simu_list_1:
    #    os.system('python fond_neutron_batch.py {}'.format(stream))
        command = 'python fond_neutron_simu_batch.py {0} {1}'.format(stream, simu)
        print(command)
        os.system(command)
    
print('Done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
