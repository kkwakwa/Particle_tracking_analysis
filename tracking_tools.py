# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 16:43:40 2018

@author: kkwakwa

This is a library of functions for importing standardised track data (from
either TrackPy or MATLABTrack), filtering it by track length, calculating the
important track statistics and plotting the data.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import trackpy as tp


def track_convert(filteredspots, dist=500):
    '''
    Takes a (hopefully standardised) Pandas dataframe of spot data and passes
    it to trackpy for spot linking. The dist option defines the maximum
    frame-to-frame distance in nm
    what is returned is the same spot dataframe, but with a new column attached
    renamed to track that includes the track ID
    '''
    tracks = tp.link(filteredspots, search_range=dist, pos_columns=['x', 'y'],
                     t_column='frame')
    tracks = tracks.rename(index=str, columns={'particle':'track_no'})
    return tracks

def read_MATLABtrack(filename):
    '''
    Takes a string that represents the path to a MATLABtrack file and returns
    a Pandas dataframe that contains the tracked spot data converted to
    the standardised header format
    '''

    trackdata = pd.read_csv(filename)
    return trackdata


def filter_tracks(tracklist, tracklength=5):
    '''
    takes a Pandas dataframe of single particle tracking data and filters
    it by minimum track length, returning the same dataframe, but with longer tracks
    '''
    tracklengths = pd.value_counts(tracklist['track_no'])
    tracklengths = tracklengths.loc[tracklengths > tracklength]
    filteredtracks = tracklist.loc[tracklist["track_no"].isin(tracklengths.index)]
    return filteredtracks


def calculate_segments(spots):
    '''
    Takes in a standardised(I hope) pandas dataframe containing information about
    the spots in track data and returns a new dataframe with per-segment track
    information. The new dataframe is fomatted:

    |track_no|start_spot|end_spot|x_displacement|y_displacement|abs_displacement|
    '''
    tracklengths = pd.value_counts(spots['track_no'])
    seglength = len(spots)-len(tracklengths)
    segs = np.zeros((seglength,6))

    counter = 0
    spotgroup = spots.groupby(['track_no'])
    for name, group in spotgroup:
        start = counter
        counter = counter + len(group) - 1
        stop = counter
        vals = group[['x','y']].values
        ids = group.index.values

        diffvals = np.diff(vals, axis=0)

        segs[start:stop,0] = name
        segs[start:stop,1] = ids[:-1]
        segs[start:stop,2] = ids[1:]
        segs[start:stop,3:5] = diffvals
        segs[start:stop,5] = np.sqrt(diffvals[:,0]**2 + diffvals[:,1]**2)

    segments = pd.DataFrame(segs, columns=['track_no',
                                                          'start_spot',
                                                          'end_spot',
                                                          'x_displacement',
                                                          'y_displacement',
                                                          'abs_displacement'])
    return segments


def MSDcalc(xypos):
    steps = len(xypos)//2
    scratch = np.zeros((steps, steps))
    for i in range(steps):
        if i ==0:
            scratch[:,i] = ((freespots[i,:] - freespots[i+1:,:])**2).sum(axis=1)
        else:
            scratch[:-i,i] = ((freespots[i,:] - freespots[i+1:,:])**2).sum(axis=1)
    freescratch = (np.ma.average(np.ma.masked_where(scratch == 0., scratch), axis=1)).data
    return freescratch


def calculate_tracks(segments):
    '''
    Takes in a standardised(I hope) pandas dataframe containing information
    about the spots in track data and returns a new dataframe with per-track
    information. The new dataframe is fomatted

    |particle|no_segments|mean_displacement|MSD|
    We are still sorting out MSD, so for now I'll skip it
    '''
    tracklengths = pd.value_counts(segments['track_no'])
    tracks = pd.DataFrame(np.zeros((len(tracklengths), 3)), columns=['track_no',
                                                            'no_segments',
                                                            'mean_displacement'])
    tracks['track_no'] = pd.value_counts(segments['track_no']).sort_index().index
    tracks['no_segments'] = tracklengths.values
    tracks['mean_displacement'] = segments.groupby(['track_no'])['abs_displacement'].mean().values
    return(tracks)


'''
Functions for plotting your track data in potentially useful ways
'''


def plot_segment_stats(segments):
    '''
    Takes the standardised Pandas dataframe containing per-segment track
    information and plots the x-displacement, y-displacement and abs
    displacement distribution while also printing out the relevant statistics
    '''
    fig = plt.figure(figsize=(10, 3))
    fig.suptitle('{} track segments'.format(len(segments)), fontsize=14)
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.hist(segments["x_displacement"], bins=100, color='b', alpha=0.5,
                                                            density=True)
    ax1.set_title('Standard Deviation = {:.2f}'.format(segments['x_displacement'].std()))
    ax1.set_xlabel('x-displacement(nm)')
    ax2 = fig.add_subplot(1,3,2)
    ax2.hist(segments["y_displacement"], bins=100, color='g', alpha=0.5, density=True)
    ax2.set_title('Standard Deviation = {:.2f}'.format(segments['y_displacement'].std()))
    ax2.set_xlabel("y-displacement(nm)")
    ax3 = fig.add_subplot(1,3,3)
    ax3.hist(segments["abs_displacement"], bins=100, color='r', alpha=0.5, density=True)
    ax3.set_title('Median = {:.2f}'.format(segments['abs_displacement'].median()))
    ax3.set_xlabel("r-displacement(nm)")
    fig.tight_layout()
    plt.subplots_adjust(top=0.8)
    #plt.show()
    return fig


def plot_track_stats(tracks):
    '''
    Takes a standardised Pandas dataframe containing per-track information
    and plots a histogram of the mean displacements
    '''
    fig = plt.figure(figsize=(7, 3))
    fig.suptitle('{} tracks'.format(len(tracks)), fontsize=14)
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.hist(tracks['mean_displacement'], bins=100, color='g', alpha=0.5, density=True)
    ax1.set_xlabel('Mean Displacement(nm)')
    ax1.set_title('Median = {:.2f}'.format(tracks['mean_displacement'].median()))
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.hist(tracks['no_segments'], bins=100, color='b', alpha=0.5, density=True)
    ax2.set_xlabel('Track lengths(no.of spots)')
    ax2.set_title('Median = {:.2f}'.format(tracks['no_segments'].median()))
    fig.tight_layout()
    plt.subplots_adjust(top=0.8)
    # plt.show()
    return fig
