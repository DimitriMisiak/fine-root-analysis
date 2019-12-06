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

for stream in stream_list:
    os.system('python fond_neutron_batch.py {}'.format(stream))
    print('Done !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')