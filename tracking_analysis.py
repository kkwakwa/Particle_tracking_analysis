# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 14:25:02 2019

@author: kkwakwa

This is a library of supporting data analysis tools that we will need for spot
tracking data analysis
"""

from scipy.optimize import curve_fit
from scipy.stats import linregress
import numpy as np


def gaussian(x, amp, cen, wid):
    return amp * np.exp(-(x-cen)**2 / wid)


def doublegaussian(x, amp1, amp2, cen1, cen2, wid1, wid2):
    return (amp1 * np.exp(-(x-cen1)**2 / wid1)) + (amp2 * np.exp(-(x-cen2)**2 / wid2))


def triplegaussian(x, amp1, amp2, amp3, cen1, cen2, cen3, wid1, wid2, wid3):
    return (amp1 * np.exp(-(x-cen1)**2 / wid1)) + (amp2 * np.exp(-(x-cen2)**2 / wid2)) + (amp3 * np.exp(-(x-cen3)**2 / wid3))


def MSD_fit(values, time=0.03):
    val_length = len(values//2)
    slope, intercept, r_value, p_value, std_err = linregress(np.log(np.arange(1, val_length+1)*0.3), np.log(values[:val_length] * 1e-6))
    return slope