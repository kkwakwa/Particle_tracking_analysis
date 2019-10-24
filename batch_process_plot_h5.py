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
import paircorrelation
import tracking_tools as trt
from glob import glob
from pandas import HDFStore
from os import path
from matplotlib.ticker import FuncFormatter
#from sys import try, except

plt.style.use('seaborn-ticks')

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


def filter_length(hf5data, tracklength=4):
    '''
    Takes an hdf5 file containing dataframes of single particle tracking data and filters
    it by minimum track length, returning the three dataframes, but with longer
    tracks
    '''
    spots = hf5data['/spots']
    segments = hf5data['/segments']
    tracks = hf5data['/tracks']
    tracklengths = pd.value_counts(spots['track_no'])
    tracklengths = tracklengths.loc[tracklengths > tracklength]
    filteredspots = spots.loc[spots["track_no"].isin(tracklengths.index)]
    filteredsegments = segments.loc[segments["track_no"].isin(tracklengths.index)]
    filteredtracks = tracks.loc[tracks["track_no"].isin(tracklengths.index)]
    return filteredspots, filteredsegments, filteredtracks


def tickformat(x):
    if int(x) == float(x):
        return str(int(x))
    else:
        return str(x)


def figurename(filename):
    tempname = path.split(filename)[-1]
    endval = tempname.find("MMStack")
    return(tempname[:endval-3])


def summary_plot(hf5filename, tracklength=4, time=0.015):
    hf5file = HDFStore(hf5filename)
    spots = hf5file['/spots']
#    segments = hf5file['/segments']
#    tracks = hf5file['/tracks']
    
    headername = figurename(hf5filename)

    
    #spots1, segments1, tracks1 = filter_length(hf5file, tracklength=2)
    spots1 = trt.filter_tracks(spots, tracklength=1)
    segments1 = trt.calculate_segments(spots1)
    tracks1 = trt.calculate_tracks(segments1)
    spots2, segments2, tracks2 = filter_length(hf5file, tracklength=4)
    
    #for plot 1
    bins, edges = np.histogram(spots.groupby('track_no').min()['frame']*time, bins=100)
    bins1, edges1 = np.histogram(spots1.groupby('track_no').min()['frame']*time, bins=100)
    bins2, edges2 = np.histogram(spots2.groupby('track_no').min()['frame']*time, bins=100)
    bins = np.cumsum(bins)
    bins1 = np.cumsum(bins1)
    bins2 = np.cumsum(bins2)
    binsd = bins/24.38**2
    binsd1 = bins1/24.38**2
    binsd2 = bins2/24.38**2
    
    spotgroup = spots.groupby('track_no')  
    tracklengths = pd.value_counts(spots['track_no'])
    segs = np.zeros((len(tracklengths), 2))
    counter = 0
    for name, group in spotgroup:
        segs[counter,:] = group[['x','y']].iloc[-1,:].values
        counter += 1
    
    spotgroup1 = spots1.groupby('track_no')
    tracklengths1 = pd.value_counts(spots1['track_no'])
    segs1 = np.zeros((len(tracklengths1), 2))
    counter1 = 0
    for name, group in spotgroup1:
        segs1[counter1,:] = group[['x','y']].iloc[-1,:].values
        counter1 += 1
    
    spotgroup2 = spots2.groupby('track_no')
    tracklengths2 = pd.value_counts(spots2['track_no'])
    segs2 = np.zeros((len(tracklengths2), 2))
    counter2 = 0
    for name, group in spotgroup2:
        segs2[counter2,:] = group[['x','y']].iloc[-1,:].values
        counter2 += 1
    
    #Calculate spot insertion time
    inserttvals = tracks1['insertion_time'].values
    inserttvals = inserttvals[~np.isnan(inserttvals)]
    insertion_ratio = tracks2['insertion_time'].isna().sum()/len(tracks2['insertion_time'])
    
#    Calculate g(r)
    segsfilter = segs2[(segs2[:,0] < 6000.) & (segs2[:,1] < 6000.)]
    try:
        g_r_2, r_2, reference_indices_2 = paircorrelation.pairCorrelationFunction_2D(segsfilter[:,0], segsfilter[:,1], 6000, 400, 20)
    except (RuntimeError):
        r_2 = np.linspace(0,400,20)
        g_r_2 = np.ones(len(r_2))
    
    #Figure parameters
    fmt = FuncFormatter(lambda x, pos: tickformat(x / 1e3))
#    fmt1 = FuncFormatter(lambda x, pos: tickformat(x / 1e6))
    fig = plt.figure(figsize=(10, 12))
    fig.suptitle(headername)

    #Number of detected plots over time at tracks = 2 spots and tracks = 5 spots
    ax1 = fig.add_subplot(4,2,1)
    ax1.plot(edges[1:], bins, 'm-', alpha = 0.5, label='min. track length = 1')
    ax1.plot(edges1[1:], bins1, 'c-', alpha = 0.5, label='min. track length = 2')
    ax1.plot(edges2[1:], bins2, 'b-', alpha = 0.5, label='min. track length = 5')
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Count ($x10^{3}$)')
#    ax1.ticklabel_format(axis='y', style='sci', scilimits=(3,3))
    ax1.yaxis.set_major_formatter(fmt)
    ax1_r = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax1_r.set_ylabel('Count per $\mu m^{2}$')  # we already handled the x-label with ax1
    ax1_r.plot(edges[1:], binsd, 'm-', alpha = 0.5, label='min. track length = 1')
    ax1_r.plot(edges1[1:], binsd1, 'c-', alpha = 0.5, label='min. track length = 2')
    ax1_r.plot(edges2[1:], binsd2, 'b-', alpha = 0.5, label='min. track length = 5')
    ax1_r.tick_params(axis='y')
    ax1_r.set_ylim(bottom=0)
    ax1.legend()
    ax1.set_title('(a)', loc='left', fontweight='bold')
    

    ax2 = fig.add_subplot(4,2,2)
    ax2.plot(r_2, g_r_2, 'b-')
    ax2.set_xlim(left=0)
    ax2.set_ylim(bottom=0)
    ax2.set_xlabel("distance (nm)")
    ax2.set_ylabel("g(r)")
    ax2.set_title("(b)", loc='left', fontweight='bold')

    
    #Mean displacement plots at track = 2 spots and track = 5 spots
    ax3 = fig.add_subplot(4,2,3)
    ax3bins = np.linspace(0,400,101)
    ax3.hist(tracks1['mean_displacement'], bins=ax3bins, color='c', alpha=0.5, label='min. track length = 2')
    ax3.hist(tracks2['mean_displacement'], bins=ax3bins, color='b', alpha=0.5, label='min. track length = 5')
    ax3.legend()
    ax3.set_xlim(left=0)
    ax3.set_ylim(bottom=0)
    ax3.set_xlabel('Mean displacement per frame (nm)')
    ax3.set_ylabel('Count ($x10^{3}$)')
#    ax3.ticklabel_format(axis='y', style='sci', scilimits=(3,3))
    ax3.yaxis.set_major_formatter(fmt)
    ax3text = ('Median = {:.2f}'.format(tracks2['mean_displacement'].median()))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax3.text(0.55, 0.55, ax3text, transform=ax3.transAxes, fontsize=12, color='b',
        verticalalignment='top', bbox=props)
    ax3.set_title('(c)', loc='left', fontweight='bold')
    
    
    #Histogram of diffusion constant
    ax4 = fig.add_subplot(4, 2, 4)
    ax4bins = np.linspace(-0.5,1.0,101)
    ax4.hist(tracks2['MSD']/4, bins=ax4bins, color='b', alpha=0.6)
    ax4.set_xlabel('Diffusion constant ($\mu m^2$/s)')
    ax4.set_ylabel('Count ($x10^{3}$)')
#    ax4.ticklabel_format(axis='y', style='sci', scilimits=(3,3))
    ax4.yaxis.set_major_formatter(fmt)
    ax4text = ('Median = {:.2f}'.format(tracks2['MSD'].median()/4))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax4.text(0.35, 0.95, ax4text, transform=ax4.transAxes, fontsize=12, color='b',
        verticalalignment='top', bbox=props)
    ax4.set_title('(d)', loc='left', fontweight='bold')
    
    
    #Histogram of track insertion times
    ax5 = fig.add_subplot(4,2,5)
    ax5bins = np.linspace(0.,2.,51)
    ax5.hist(inserttvals, bins=ax5bins, color='b', alpha=0.5)
    ax5.set_yscale('log')
    ax5.set_xlim(left=0)
    ax5.set_ylim(bottom=0)
    ax5.set_xlabel('Time to insertion (s)')
    ax5.set_ylabel('Count')
    ax5.minorticks_off()
    #ax5.ticklabel_format(axis='y', style='sci', scilimits=(3,3))
    ax5text = ('{:.1f}% tracks inserted \n{:.3f}s median insertion time'.format((1-insertion_ratio)*100, tracks2['insertion_time'].median()))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax5.text(0.25, 0.95, ax5text, transform=ax5.transAxes, fontsize=12, color='b',
        verticalalignment='top', bbox=props)
    ax5.set_title('(e)', loc='left', fontweight='bold')
    
    
    #Histogram of track duration in seconds
    ax6 = fig.add_subplot(4,2,6)
    ax6bins = np.linspace(0,5,51)
    ax6.hist(tracks2['no_segments']*.015, bins=ax6bins, color='b', alpha=0.5)
    ax6.set_yscale('log')
    ax6.set_xlim(left=0)
    #ax6.set_ylim(bottom=0)
    #ax6.ticklabel_format(axis='y', style='sci', scilimits=(3,3))
    ax6.set_xlabel('Track lengths (s)')
    ax6.set_ylabel('Count')
    ax6.minorticks_off()
    ax6text=('Median = {:.2f}s'.format(tracks2['no_segments'].median()*time))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax6.text(0.55, 0.95, ax6text, transform=ax6.transAxes, fontsize=12, color='b',
        verticalalignment='top', bbox=props)
    ax6.set_title('(f)', loc='left', fontweight='bold')
    
    
    #Histogram of mean track intensity
    ax7 = fig.add_subplot(4,2,7)
    ax7bins = np.linspace(0,1000,101)
    ax7.hist(spotgroup.mean()['intensity'], bins=ax7bins, color='m', alpha=0.5, label='min. track length = 1')
    ax7.hist(spotgroup1.mean()['intensity'], bins=ax7bins, color='c', alpha=0.5, label='min. track length = 2')
    ax7.hist(spotgroup2.mean()['intensity'], bins=ax7bins, color='b', alpha=0.5, label='min. track length = 5')
    ax7.set_xlim(left=0)
    ax7.set_ylim(bottom=0)
    ax7.set_xlabel('Mean track intensity(photons)')
    ax7.set_ylabel('Count ($x10^{3}$)')
    ax7.legend()
#    ax7.ticklabel_format(axis='y', style='sci', scilimits=(3,3))
    ax7.yaxis.set_major_formatter(fmt)
    ax7text = ('Median = {:.2f}'.format(spotgroup2.mean()['intensity'].median()))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax7.text(0.55, 0.55, ax7text, transform=ax7.transAxes, fontsize=12, color='b',
        verticalalignment='top', bbox=props)
    ax7.set_title('(g)', loc='left', fontweight='bold')
    
    
    #2D histogram of final track positions
    ax8 = fig.add_subplot(4,2,8)
    plt.hist2d(segs2[:,0], segs2[:,1], bins=(500,500), cmap='Greys')
    ax8.set_aspect(aspect=1)
    ax8.set_yticklabels([])
    ax8.set_xticklabels([])
    ax8.set_title('(h)', loc='left', fontweight='bold')
    plt.colorbar()
    
    
    fig.tight_layout()
    plt.subplots_adjust(top=0.92)
    return fig
    #fig.show()


directory = 'C:/Data/Fluorescence/Me/2019_10_10/'
search = '**/*trackpy.h5'
#search = '*trackpy.csv'

#directory = 'C:/Data/Fluorescence/Me/'
#search = '**/**/*trackpy.h5'

#search = '*all_ts.csv'

filelist = glob(directory+search)


for file in filelist:
    if "bleach" in file:
        filelist.remove(file)
        
for file in filelist:
    if "bleach" in file:
        filelist.remove(file)
        
print(filelist)

for file in filelist:
    print(file)
    
    if "fast" in file:
        ttime = 0.005
    elif "slow" in file:
        ttime = 0.030
    else:
        ttime= 0.015
    
    longplotfile = (file[:-6]+'summary_plot.png')
    longplotfile1 = (file[:-6]+'summary_plot.svg')

#    spotdata = HDFStore(file)
#    print('loaded')

    longfig = summary_plot(file, time=ttime)


    print('Saving everything........')
    longfig.savefig(longplotfile, dpi=300, frameon=True)
    longfig.savefig(longplotfile1)
    print('Done!')
    plt.close(longfig)
print('All Done!')
