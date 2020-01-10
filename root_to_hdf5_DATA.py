#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 10:09:10 2020

@author: misiak
"""

import os
import uproot
import pandas as pd
from tqdm import tqdm

stream_config = {
    "Background": [
    	"tg18l005",
    	"tg27l000",
    	"tg28l000"
    ],
    "Calibration": [
    	"tg17l007",
    	"tg19l010",
    	"tg20l000",
    	"tg21l000"
    ]
}

detector = "RED80"
data_dir = '/home/misiak/Data/data_run57_neutron/Data'
simu_dir = '/home/misiak/Data/data_run57_neutron/SimuCoinc'


# detector, mode, stream, file path
stream_info_list = list()

for mode, stream_list in stream_config.items():
    
    for stream in stream_list:
        
        dir_path = '/'.join([data_dir, stream, detector])
        
        for fname in os.listdir(dir_path):
            
            if 'ProcessedData' in fname:
                
                file_path = '/'.join([dir_path, fname])
                stream_info = [
                    detector,
                    mode,
                    stream,
                    file_path
                ]
                stream_info_list.append(stream_info)


# def decode_byte_column_name(df):
#     """
#     Decode the columns name of bytes type into str type.
#     Parameters
#     ----------
#     df : pandas.DataFrame

#     Returns
#     -------
#     None.
#     """
#     decode_dict = {b:str(b, "utf-8") for b in df.columns if isinstance(b, bytes)}
#     df.rename(columns = decode_dict, inplace=True)    
#     return


def root_to_df(root):

    # joining trees for triggering events
    df_raw = root['EventTree_trig_Normal_raw'].pandas.df()
    df_filt = root['EventTree_trig_Normal_filt'].pandas.df()
    df_filt_decor = root['EventTree_trig_Normal_filt_decor'].pandas.df()
    
    df_list = [df_raw, df_filt, df_filt_decor]
    suffix_list = ['_raw', '_filt', '_filt_decor']
    for df, suffix in zip(df_list, suffix_list):
        # decode_byte_column_name(df)
        df.rename(
            columns = {name:name+suffix for name in df.columns},
            inplace=True
        )

    
    df = df_raw.join([df_filt, df_filt_decor])
    
    
    # new columns for general info in the RunTree_Normal tree
    run_tree = root['RunTree_Normal']
    df['frequency'] = run_tree['f_max_heat'].array()[0]
    df['time_window'] = run_tree['TimeWindow_Heat'].array()[0]
    df['maintenance_cycle'] = run_tree['MaintenanceCycle'].array()[0]
    df['maintenance_duration'] = run_tree['MaintenanceDuration'].array()[0]
    
    polarization = run_tree['Polar_Ion'].array()[0]
    for polar, label in zip(polarization, ['A', 'B', 'C', 'D']):
        df['polar_'+label] = polar
    
    return df


# df_list = list()

save_dir = '/home/misiak/Analysis/neutron_background'
save_path = '/'.join([save_dir, 'data.h5'])


with pd.HDFStore(save_path, mode='w') as store:
    
    chunk_start = 0
    for detector, mode, stream, file_path in tqdm(stream_info_list):
        
        # opening root file
        root = uproot.open(file_path)
        
        # creating a DataFrame from the info in the RootDirectory
        df = root_to_df(root)
        
        # converting bytes to str
        str_df = df.select_dtypes('O')
        str_df = str_df.stack().str.decode('utf-8').unstack()
        for col in str_df:
            df[col] = str_df[col]
            
        # creating columns for more general info
        df['detector'] = detector
        df['source'] = mode
        df['stream'] = stream   
        
        # correct indexing
        chunk_end = chunk_start + df.shape[0]
        df.index = range(chunk_start, chunk_end)
        chunk_start = chunk_end
        
        # df_list.append(df)
        store.append('data', df, min_itemsize=11)
