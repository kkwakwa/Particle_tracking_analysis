# -*- coding: utf-8 -*-
"""
@author: kkwakwa

This is a library of functions extracting bleach data from image stacks, then
running analysis on those to get both single molecule bleach times and count
the steps in multi-molecule bleaching
"""


import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from skimage import io
from skimage.feature import peak_local_max
from skimage.filters import threshold_otsu


def steps_extract(filename):
    bleach_stack = io.imread(filename)
    max_intensity = bleach_stack.max(0)
    thresh = threshold_otsu(max_intensity)
    peaks = peak_local_max(max_intensity, min_distance=4, threshold_abs=thresh)
    traces = []
    no = 0
    for spot in peaks:
        xmin = spot[0] - 2
        xmax = spot[0] + 3
        ymin = spot[1] - 2
        ymax = spot[1] + 3
        trace = bleach_stack[:, xmin:xmax, ymin:ymax].sum(axis=(1, 2))
        time = np.arange(len(trace))*0.1
        traces.append(pd.Series(trace, index=time, name="track_"+repr(no)))
        no += 1
    peaklist = pd.concat(traces, axis=1)
    outfile = filename[:-4]+'-out.csv'
    peaklist.to_csv(path_or_buf=outfile)
