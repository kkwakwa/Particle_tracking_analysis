# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:16:37 2018

@author: kkwakwa
"""

import spot_tools as spt
from glob import glob


directory = 'D:/Data/Fluorescence/Kwasi/2018_12_21/'
search = '**/*_ts.csv'


filelist = glob(directory+search)
print(filelist)


for file in filelist:
    print('Processing: '+file)
    
    
    spottimesplot = (file[:-23]+'ts_spottimes.png')
    spotstatsplot = (file[:-23]+'ts_spotstats.png')
    
    spotfile = spt.readThunderStorm(file)
    filteredspots = spt.filter_spots(spotfile)
    
    
    spottimesfig = spt.plotspottimes(filteredspots)
    spotsstatsfig = spt.plotspotstats(filteredspots)
    
    
    print('Saving everything........')
    spottimesfig.savefig(spottimesplot, dpi=300, frameon=True)
    spotsstatsfig.savefig(spotstatsplot, dpi=300, frameon=True)
    print('Done!')