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
from scipy.stats import linregress


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
    tracks = tracks.rename(index=str, columns={'particle': 'track_no'})
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
    Takes a Pandas dataframe of single particle tracking data and filters
    it by minimum track length, returning the same dataframe, but with longer
    tracks
    '''
    tracklengths = pd.value_counts(tracklist['track_no'])
    tracklengths = tracklengths.loc[tracklengths > tracklength]
    filteredtracks = tracklist.loc[tracklist["track_no"].isin(tracklengths.index)]
    return filteredtracks


def MSDcalc(xypos):
    '''
    Takes a two column numpy array where the first column is all particle
    x-positions and the second column is particle y-positions, returns a single
    column one step shorter than the column length of the input column. Each
    position in that column is the MSD step for that length
    '''
    steps = len(xypos) - 1
    scratch = np.zeros((steps, steps))
    for i in range(steps):
        if i == 0:
            scratch[:, i] = ((xypos[i, :] - xypos[i+1:, :])**2).sum(axis=1)
        else:
            scratch[:-i, i] = ((xypos[i, :] - xypos[i+1:, :])**2).sum(axis=1)
    freescratch = (np.ma.average(np.ma.masked_where(scratch == 0., scratch), axis=1)).data
    return freescratch


def MSD_fit(values, time=0.03):
    val_length = len(values//2)
    slope, intercept, r_value, p_value, std_err = linregress(np.log(np.arange(1, val_length+1)*time), np.log(values[:val_length]))
    return slope


def odist(xypos):
    '''
    takes a two column numpy array where the first column is particle
    x-positions and the second column is particle y-positions.

    Returns a single column that is one step shorter than the input column.
    Each position is a measure of the square of the length of displacement from
    the starting position
    '''

    return((xypos[-1, :] - xypos[:-1, :])**2).sum(axis=1)


def calculate_segments(spots):
    '''
    Takes in a standardised(I hope) pandas dataframe containing information
    about the spots in track data and returns a new dataframe with per-segment
    track information. The new dataframe is fomatted:

    |track_no|start_frame|start_spot|end_spot|x_displacement|y_displacement|
    abs_displacement|MSDs|origin_dist|
    '''
    tracklengths = pd.value_counts(spots['track_no'])
    seglength = len(spots)-len(tracklengths)
    segs = np.zeros((seglength, 9))

    counter = 0
    spotgroup = spots.groupby(['track_no'])
    for name, group in spotgroup:
        start = counter
        counter = counter + len(group) - 1
        stop = counter
        vals = group[['x', 'y']].values
        ids = group.index.values

        diffvals = np.diff(vals, axis=0)

        segs[start:stop, 0] = name
        segs[start:stop, 1] = group['frame'].values[:-1]
        segs[start:stop, 2] = ids[:-1]
        segs[start:stop, 3] = ids[1:]
        segs[start:stop, 4:6] = diffvals
        segs[start:stop, 6] = np.sqrt(diffvals[:, 0]**2 + diffvals[:, 1]**2)
        segs[start:stop, 7] = MSDcalc(vals)
        segs[start:stop, 8] = odist(vals)

    segments = pd.DataFrame(segs, columns=['track_no',
                                           'start_frame',
                                           'start_spot',
                                           'end_spot',
                                           'x_displacement',
                                           'y_displacement',
                                           'abs_displacement',
                                           'MSDs',
                                           'origin_dist'])
    return segments


def calculate_tracks(segments):
    '''
    Takes in a standardised(I hope) pandas dataframe containing information
    about the spots in track data and returns a new dataframe with per-track
    information. The new dataframe is fomatted

    |particle|start_frame|no_segments|mean_displacement|MSD|
    We are still sorting out MSD, so for now I'll skip it
    '''
    tracklengths = pd.value_counts(segments['track_no'])
    counter = 0
    tracks = np.zeros((len(tracklengths), 6))
    seggroup = segments.groupby('track_no')

    for name, group in seggroup:
        tracks[counter, 0] = name
        tracks[counter, 1] = len(group)
        tracks[counter, 2] = group['start_frame'].min()
        tracks[counter, 3] = group['abs_displacement'].mean()
        tracks[counter, 4] = MSD_fit(group['MSDs'].values, time=0.015)
        temps = np.where(np.sqrt(group['origin_dist'].values) <= 70)
        if len(temps[0]) == 0:
            tracks[counter, 5] = np.NaN
        else:
            tracks[counter, 5] = ((temps[0][0] +1) *0.015 )
        counter += 1


    finaltracks = pd.DataFrame(tracks, columns=['track_no',
                                                'no_segments',
                                                'start_frame',
                                                'mean_displacement',
                                                'MSD',
                                                'insertion_time'])
    return(finaltracks)


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
    ax2 = fig.add_subplot(1, 3, 2)
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
    fig = plt.figure(figsize=(10, 3))
    fig.suptitle('{} tracks'.format(len(tracks)), fontsize=14)
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.hist(tracks['mean_displacement'], bins=100, color='g', alpha=0.5, density=True)
    ax1.set_xlabel('Mean Displacement(nm)')
    ax1.set_title('Median = {:.2f}'.format(tracks['mean_displacement'].median()))
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.hist(tracks['no_segments'], bins=100, color='b', alpha=0.5, density=True)
    ax2.set_xlabel('Track lengths(no.of spots)')
    ax2.set_title('Median = {:.2f}'.format(tracks['no_segments'].median()))
    ax3 = fig.add_subplot(1, 3, 3)
    ax3.hist(tracks['MSD'], bins=100, color='g', alpha=0.6, density=True)
    ax3.set_xlabel('MSD fit($\mu m^2$/sec)')
    fig.tight_layout()
    plt.subplots_adjust(top=0.8)
    # plt.show()
    return fig


def plot_msd_dist(tracks):
    '''
    Takes a standardised Pandas dataframe containing per-track information
    and plots a histogram of the MSD
    '''
    fig = plt.figure(figsize=(10, 6))
    plt.hist(tracks['MSD'], bins=100, color='g', alpha=0.6, density=True)
    plt.xlabel('MSD fit($\mu m^2$/sec)')
    plt.tight_layout()
    return fig


def plot_msds(tracks, time=0.015):
    '''
    Takes a standardised Pandas dataframe containing per-track information and
    plots each individual MSD
    '''
    drawlines = []
    tracksgroup = tracks.groupby("track_no")

    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_subplot(1, 1, 1)
    for name, group in tracksgroup:
        drawlines.append(np.arange(1, len(group["MSDs"].values)+1)*time)
        drawlines.append(group["MSDs"].values * 1e-6)
    ax1.loglog(*drawlines, color='0.2', linewidth=0.5, alpha=0.1)
    ax1.set_xlabel('time(s)')
    ax1.set_ylabel('Mean Squared Displacement($\mu$m$^{2}$)')
    ax1.set_title('labelled MWT')
    fig.tight_layout()
    return fig


def plot_final_positions(tracks, time=0.015):
    '''
    Takes a standardised Pandas dataframe containing per-track information and
    plots distances from the final position for each track, which is a good way
    to tell how quickly spots are settling down
    '''
    drawlines = []
    tracksgroup = tracks.groupby("track_no")
    for name, group in tracksgroup:
        drawlines.append(np.arange(1, len(group["origin_dist"].values)+1)*time)
        drawlines.append(np.sqrt(group["origin_dist"].values))

    fig = plt.figure(figsize=(10, 3))
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(*drawlines, color='0.6', linewidth=0.5, alpha=0.5)
    ax1.set_xlabel('time since appearance(s)')
    ax1.set_ylabel('Distance from final step($nm$)')
    ax1.set_title("Labelled GFP-TMH minus DTT")
    fig.tight_layout()
    return fig
