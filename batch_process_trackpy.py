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

directory = 'C:/Data/Fluorescence/Me/2019_08_14/'
search = '**/*ts.csv'
#search = '*ts.csv'

#search = '*all_ts.csv'

filelist = glob(directory+search)
print(filelist)


for file in filelist:
    print(file)

    trackfile = (file[:-6]+'trackpy.csv')
    segdatafile = (file[:-6]+'trackpy_segments.csv')
    trackdatafile = (file[:-6]+'trackpy_tracks.csv')
    spotdatafile = (file[:-6]+'trackpy_spotcount.png')
    spotplotfile = (file[:-6]+'trackpy_spotstats.png')
    spotposplotfile = (file[:-6]+'trackpy_spotpos.png')
    segplotfile = (file[:-6]+'trackpy_segmentstats.png')
    trackplotfile = (file[:-6]+'trackpy_trackstats.png')
    insertiondistfile = (file[:-6]+'trackpy_insertiondist.png')
    disttimeplotfile = (file[:-6]+'trackpy_dist_time.png')

    spotdata = spt.readThunderStorm(file)
    print('loaded')
    spotdata = spt.filter_spots(spotdata)
    spotcount = spt.plotspottimes(spotdata)

    trackdata = trt.track_convert(spotdata)
    trackdata1 = trt.filter_tracks(trackdata, tracklength=3)
    #trackdata1 = trackdata

    spotfig = spt.plotspotstats(trackdata1)
    spotdist = spt.plotspotpositions(trackdata1)

    segmentstats = trt.calculate_segments(trackdata1)
    segmentfig = trt.plot_segment_stats(segmentstats)


    trackstats = trt.calculate_tracks(segmentstats)
    trackfig = trt.plot_track_stats(trackstats)
    positiondistfig = trt.plot_final_position_dist(trackstats)
    finaldisplacementdistfig = trt.plot_displacement_time(trackstats)
    

    print('Saving everything........')
    trackdata1.to_csv(trackfile)
    spotcount.savefig(spotdatafile, dpi=300, frameon=True)
    spotfig.savefig(spotplotfile, dpi=300, frameon=True)
    spotdist.savefig(spotposplotfile, dpi=300, frameon=True)
    segmentfig.savefig(segplotfile, dpi=300, frameon=True)
#    segmentstats.to_csv(segdatafile)
    trackfig.savefig(trackplotfile, dpi=300, frameon=True)
    positiondistfig.savefig(disttimeplotfile, dpi=300, frameon=True)
    finaldisplacementdistfig.savefig(insertiondistfile, dpi=300, frameon=True)
    print('Done!')

print('All Done!')
