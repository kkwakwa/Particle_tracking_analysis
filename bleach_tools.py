# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 16:43:40 2018

@author: kkwakwa
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def calculate_segments(spots):
    '''
    Takes in a standardised(I hope) pandas dataframe containing information about
    the spots in track data and returns a new dataframe with per-segment track
    information. The new dataframe is fomatted:
    
    |track|start_spot|end_spot|x_displacement|y_displacement|abs_displacement|
    '''
    tracklengths = pd.value_counts(spots['track_no'])
    seglength = len(spots)-len(tracklengths)
    segments = pd.DataFrame(np.zeros((seglength,6)), columns=['track','start_spot',
                                                              'end_spot','x_displacement',
                                                              'y_displacement','abs_displacement'])
    counter = 0
    spotgroup = spots.groupby(['track_no'])
    for name, group in spotgroup:
        start = counter
        counter = counter + len(group) - 1
        stop = counter
        vals = group[['x','y']].values
        ids = group.index.values

        diffvals = np.diff(vals, axis=0)

        segments.iloc[start:stop,0] = name
        segments.iloc[start:stop,1] = ids[:-1]
        segments.iloc[start:stop,2] = ids[1:]
        segments.iloc[start:stop,3:5] = diffvals
        segments.iloc[start:stop,5] = np.sqrt(diffvals[:,0]**2 + diffvals[:,1]**2)
    return segments