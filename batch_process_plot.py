# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 12:01:14 2019

@author: kkwakwa
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 13 10:23:58 2019

@author: kkwakwa
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import spot_tools as spt
import tracking_tools as trt
from glob import glob

def summary_plot(spots, tracklength=4):
    segments = trt.calculate_segments(spots)
    tracks = trt.calculate_tracks(segments)
    
    spots1 = trt.filter_tracks(spots, tracklength=tracklength)
    segments1 = trt.calculate_segments(spots1)
    tracks1 = trt.calculate_tracks(segments1)
    
    tracklengths = pd.value_counts(spots['track_no'])
    segs = np.zeros((len(tracklengths), 2))
    counter = 0
    spotgroup = spots.groupby('track_no')
    for name, group in spotgroup:
        segs[counter,:] = group[['x','y']].iloc[-1,:].values
        counter += 1
    
    inserttvals = tracks1['insertion_time'].values
    inserttvals = inserttvals[~np.isnan(inserttvals)]
    insertion_ratio = tracks1['insertion_time'].isna().sum()/len(tracks['insertion_time'])
    
    fig = plt.figure(figsize=(10, 12))
    #fig.suptitle('15ms per frame for this experiment', fontsize=16)

    #Mean displacement plots at track = 2 spots and track = 5 spots
    ax1 = fig.add_subplot(3,2,1)
    ax1.hist(tracks['mean_displacement'], bins=100, color='g', alpha=0.5, label='minimum track length = 2')
    ax1.hist(tracks1['mean_displacement'], bins=100, color='r', alpha=0.5, label='minimum track length = 5')
    ax1.legend()
    ax1.set_xlabel('Mean Displacement per frame(nm)')
    ax1.set_title('Median = {:.2f}'.format(tracks['mean_displacement'].median()))
    
    ax2 = fig.add_subplot(3, 2, 2)
    ax2.hist(tracks1['MSD']/4, bins=100, color='g', alpha=0.6)
    ax2.set_xlabel('Diffusion constant($\mu m^2$/sec)')
    ax2.set_title('Median = {:.2f}'.format(tracks1['MSD'].median()/4))
    
    ax3 = fig.add_subplot(3,2,3)
    ax3.hist(inserttvals, bins=50, color='b', alpha=0.3)
    ax3.set_xlabel('time to insertion(s)')
    ax3.set_title('{:.1f}% tracks inserted, {:.3f}s median insertion time'.format((1-insertion_ratio)*100, tracks['insertion_time'].median()))
    
    ax4 = fig.add_subplot(3,2,4)
    ax4.hist(tracks1['no_segments']*.015, bins=100, color='b', alpha=0.5)
    ax4.set_xlabel('Track lengths(s)')
    ax4.set_title('Median = {:.2f}'.format(tracks1['no_segments'].median()))
    
    ax5 = fig.add_subplot(3,2,5)
    plt.hist2d(segs[:,0], segs[:,1], bins=(500,500), cmap='Greys')
    plt.colorbar()
    
    ax6 = fig.add_subplot(3,2,6)
    ax6.hist(tracks['start_frame']*0.015, bins=100, color='y', cumulative=True)
    ax6.set_xlabel('Time(s)')
    ax6.set_ylabel('No. of tracks detected')
    
    fig.tight_layout()
    #plt.subplots_adjust(top=0.92)
    return fig


directory = 'C:/Data/Fluorescence/Me/2019_08_30/'
search = '**/*trackpy.csv'
#search = '*trackpy.csv'

#search = '*all_ts.csv'

filelist = glob(directory+search)
print(filelist)

for file in filelist:
    print(file)

    longplotfile = (file[:-6]+'summary_plot.png')
#    longplotfile = (file[:-6]+'summary_plot.svg')

    spotdata = pd.read_csv(file, index_col=0)
    print('loaded')

    longfig = summary_plot(spotdata)


    print('Saving everything........')
    longfig.savefig(longplotfile, dpi=300, frameon=True)
    print('Done!')

print('All Done!')
