# -*- coding: utf-8 -*-
"""
Created on Sun Sep  8 16:27:08 2019

@author: kkwakwa

batch processing script for taking existing trackpy files, processing step and 
track tables and then saving all of it to an HDF5 file
"""

import tracking_tools as trt
from glob import glob
import pandas as pd
from pandas import HDFStore

directory = 'C:/Data/Fluorescence/Me/2019_07_26/'
search = '**/*trackpy.csv'

filelist = glob(directory+search)
print(filelist)

for file in filelist:
    print(file)

    outfile = (file[:-4]+'.h5')
    
    if "fast" in file:
        ttime = 0.005
    elif "slow" in file:
        ttime = 0.030
    else:
        ttime= 0.015
    
    spotdata = pd.read_csv(file, index_col=0)
    #decimals = pd.Series([1, 1, 1, 1, 1], index=['x', 'y', 'intensity', 
    #                                             'sigma1', 'sigma2', ])
    decimals = pd.Series([1, 1, 1, 1, 1], index=['x', 'y', 'intensity', 
                                                 'sigma', 'uncertainty', ])
    spotdata = spotdata.round(decimals)
    spotdata1 = trt.filter_tracks(spotdata, tracklength=2)
    segdata = trt.calculate_segments(spotdata1)
    trackdata = trt.calculate_tracks(segdata, time=ttime)
    
    hdf = HDFStore(outfile)
    hdf.put('spots', spotdata, format='table', data_columns=True)
    hdf.put('segments', segdata, format='table', data_columns=True)
    hdf.put('tracks', trackdata, format='table', data_columns=True)
    hdf.close()
    
    print("File saved")
    
print("Directory done")