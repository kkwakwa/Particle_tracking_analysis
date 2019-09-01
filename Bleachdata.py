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
from skimage.filters import median

from glob import glob

directory = 'D:/Data/Fluorescence/Roxie/190606/'
search = '**/*ome.tif'
#search = '*ts.csv'

#search = '*all_ts.csv'

filelist = glob(directory+search)
print(filelist)


def steps_extract(filename):
    bleach_stack = io.imread(filename)
    max_intensity = bleach_stack[:1000].max(0)
    max_intensity_1 = median(max_intensity)
    thresh = threshold_otsu(max_intensity_1)
    peaks = peak_local_max(max_intensity_1, min_distance=4, threshold_abs=thresh+0.0001)
    traces = []
    no = 0
    for spot in peaks:
        xmin = spot[0] - 3
        xmax = spot[0] + 4
        ymin = spot[1] - 3
        ymax = spot[1] + 4
        trace = bleach_stack[:, xmin:xmax, ymin:ymax].sum(axis=(1, 2))
        time = np.arange(len(trace))*0.1
        traces.append(pd.Series(trace, index=time, name="track_"+repr(no)))
        no += 1
    peaklist = pd.concat(traces, axis=1)
    #background = peaklist.min(axis=1)
    #peaklist1 = peaklist.subtract(background, axis=0)
    outfile = filename[:-4]+'-out.csv'
    peaklist.to_csv(path_or_buf=outfile)

for file in filelist:
    print(file)
    steps_extract(file)