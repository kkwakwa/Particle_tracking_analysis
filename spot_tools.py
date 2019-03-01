# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 12:36:58 2018

@author: kkwakwa

This is a series of functions for importing spot data(currently from
ThunderSTORM), filtering out spots with based on sigma and uncertainty, and then
plotting out spot statistics

"""

#import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def readThunderStorm(filename):
    '''
    Takes a string that represents the path to a ThunderSTORM file and returns
    a Pandas dataframe that contains the most important parts of the dataframe
    using a standardised header format
    '''
    spotdata = pd.read_csv(filename, index_col=0, usecols=['id', 'frame', 'x [nm]',
                                                           'y [nm]', 'sigma [nm]',
                                                           'intensity [photon]',
                                                           'uncertainty_xy [nm]'])
    spotdata = spotdata.rename(index=str, columns={'x [nm]':'x', 'y [nm]':'y',
                                                   'sigma [nm]': 'sigma',
                                                   'intensity [photon]':'intensity',
                                                   'uncertainty_xy [nm]':'uncertainty'})
    return spotdata


def filter_spots(spotdata, sigma_min=50, sigma_max=250, unc_min=5, unc_max=65):
    '''
    Takes a (hopefully standardised) Pandas dataframe of spot data and returns
    it filtered by sigma and uncertainty.
    '''
    filteredspots = spotdata.loc[(spotdata['sigma'] > sigma_min) &
                                 (spotdata['sigma'] < sigma_max) &
                                 (spotdata['uncertainty'] > unc_min) &
                                 (spotdata['uncertainty'] < unc_max)]
    return(filteredspots)


'''
Functions for plotting your spot data in potentially useful ways
'''


def plotspottimes(spotdata):
    '''
    Takes a (hopefully standardised) Pandas dataframe of spot data and returns
    a histogram of spot count against frame number, so we can see if anything
    changes over time
    '''
    f = plt.figure(figsize=(10, 3))
    ax1 = f.add_subplot(1, 1, 1)
    ax1.hist(spotdata['frame'], bins=100, color='y')
    ax1.set_xlabel('Frame')
    ax1.set_ylabel('No. of spots detected')
    f.tight_layout()
    #plt.show()
    return f


def plotspotstats(spotdata):
    '''
    Takes a (hopefully standardised) Pandas dataframe of spot data and plots
    the distributions of sigma, uncertainty and intensity. The plots also
    include the median value of each of those.
    '''
    f = plt.figure(figsize=(10, 3))
    f.suptitle('{} spots'.format(len(spotdata)), fontsize=14)
    ax1 = f.add_subplot(1, 3, 1)
    ax1.set_title('median = {:.2f}'.format(spotdata['sigma'].median()))
    ax1.hist(spotdata['sigma'], bins=100, color='b', alpha=0.6, density=True)
    plt.xlabel('Sigma(nm)')
    ax2 = f.add_subplot(1, 3, 2)
    ax2.set_title('median = {:.2f}'.format(spotdata['uncertainty'].median()))
    ax2.hist(spotdata['uncertainty'], bins=100, color='y', alpha=0.6, density=True)
    plt.xlabel('Uncertainty(nm)')
    ax3 = f.add_subplot(1, 3, 3)
    ax3.set_title('median = {:.2f}'.format(spotdata['intensity'].median()))
    ax3.hist(spotdata['intensity'], bins=100, color='g', alpha=0.6, density=True)
    plt.xlabel('Intensity(Photons)')
    f.tight_layout()
    plt.subplots_adjust(top=0.8)
    #plt.show()
    return f
