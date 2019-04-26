# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 10:23:58 2019

@author: kkwakwa
"""

#import numpy as np
#import matplotlib.pyplot as plt
#import pandas as pd
import spot_tools as spt
import tracking_tools as trt
from glob import glob

directory = 'D:/Data/Fluorescence/Kwasi/2019_04_16/'
search = '**/*_ts.csv'


filelist = glob(directory+search)
print(filelist)


for file in filelist:
    print(file)
    
    trackfile = (file[:-6]+'trackpy.csv')
    segdatafile = (file[:-6]+'trackpy_segments.csv')
    trackdatafile = (file[:-6]+'trackpy_tracks.csv')
    spotdatafile = (file[:-6]+'trackpy_spotcount.png')
    spotplotfile = (file[:-6]+'trackpy_spotstats.png')
    segplotfile = (file[:-6]+'trackpy_segmentstats.png')
    trackplotfile = (file[:-6]+'trackpy_trackstats.png')
    
    spotdata = spt.readThunderStorm(file)
    spotdata = spt.filter_spots(spotdata)
    spotcount = spt.plotspottimes(spotdata)
    
    trackdata = trt.track_convert(spotdata)
    trackdata1 = trt.filter_tracks(trackdata, tracklength=2)
    #trackdata1 = trackdata
    
    spotfig = spt.plotspotstats(trackdata1)
    
    
    segmentstats = trt.calculate_segments(trackdata1)
    segmentfig = trt.plot_segment_stats(segmentstats)
    
    
    trackstats = trt.calculate_tracks(segmentstats)
    trackfig = trt.plot_track_stats(trackstats)
    
    print('Saving everything........')
    trackdata1.to_csv(trackfile)
    spotcount.savefig(spotdatafile, dpi=300, frameon=True)
    spotfig.savefig(spotplotfile, dpi=300, frameon=True)
    segmentfig.savefig(segplotfile, dpi=300, frameon=True)
#    segmentstats.to_csv(segdatafile)
    trackfig.savefig(trackplotfile, dpi=300, frameon=True)
    print('Done!')
    
print('All Done!')