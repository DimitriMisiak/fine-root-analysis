#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 10:59:21 2019

@author: misiak
"""

import uproot
import numpy as np

A = uproot.open('/home/misiak/Data/data_run57_neutron/SimuCoinc/tg20l000/RED80/Flat_0.0000_2.0000_ER_5/SIMUCOINC_ProcessedData_tg20l000_S00_RED80_Flat_0.0000_2.0000_ER_5_ChanTrig0.root')


B = A['tsimucoinc']

C = A['tree_Input_SimuOnly_Energy']