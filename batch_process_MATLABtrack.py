# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 14:17:19 2018

@author: kkwakwa

This is a script for batch processing tracking data from MATLABtrack data
The idea is to process the data and save out a series of graphs with the results
that I'm currently interested in

"""

#import numpy as np
#import matplotlib.pyplot as plt
#import pandas as pd
#import spot_tools as spt
import tracking_tools as trt
from glob import glob

directory = 'D:/Data/Fluorescence/Kwasi/2018_12_21/'
search = '**/*_matlabtrack.csv'


filelist = glob(directory+search)
print(filelist)


for file in filelist:
    print(file)
        
    segdatafile = (file[:-4]+'_trackpy_segments.csv')
    trackdatafile = (file[:-4]+'_trackpy_tracks.csv')
    spotdatafile = (file[:-4]+'_trackpy_spotcount.png')
    spotplotfile = (file[:-4]+'_trackpy_spotstats.png')
    segplotfile = (file[:-4]+'_trackpy_segmentstats.png')
    trackplotfile = (file[:-4]+'_trackpy_trackstats.png')
        
 
    trackdata1 = trt.read_MATLABtrack(file)
    
                
    segmentstats = trt.calculate_segments(trackdata1)
    segmentfig = trt.plot_segment_stats(segmentstats)
        
        
    trackstats = trt.calculate_tracks(segmentstats)
    trackfig = trt.plot_track_stats(trackstats)
        
    
    print('Saving everything........')
    segmentfig.savefig(segplotfile, dpi=300, frameon=True)
    segmentstats.to_csv(segdatafile)
    trackfig.savefig(trackplotfile, dpi=300, frameon=True)
    print('Done!')
        
print('All Done!')
