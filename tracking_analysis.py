# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:25:02 2019

@author: kkwakwa

This is a library of supporting data analysis tools that we will need for spot
tracking data analysis
"""

import numpy as np
from scipy.optimize import curve_fit


def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)


def doublegaussian(x, amp1, amp2, cen1, cen2, wid1, wid2):
    return (amp1 * np.exp(-(x-cen1)**2 / wid1)) + (amp2 * np.exp(-(x-cen2)**2 / wid2))


def triplegaussian(x, amp1, amp2, amp3, cen1, cen2, cen3, wid1, wid2, wid3):
    return (amp1 * np.exp(-(x-cen1)**2 / wid1)) + (amp2 * np.exp(-(x-cen2)**2 / wid2)) + (amp3 * np.exp(-(x-cen3)**2 / wid3))


def gaussian_fit_data(data):
    '''
        takes a single column numpy array representing a range of values and
        then attempts to fit a double gaussian

        Returns an array consisting of [amp1, amp2, cen1, cen2, wid1, wid2]

        standard deviation = sqrt(wid/2)
    '''
    hist, edges = np.histogram(data, bins=100)
    edges1 = edges[:-1] + np.diff(edges)[0]
    init_vals = [.5, .5, 50, 250, 10, 10]  # for [amp1,amp2, cen1, cen2, wid1, wid2]
    best_vals, covar = curve_fit(doublegaussian, edges1, hist, p0=init_vals)
    return best_vals
