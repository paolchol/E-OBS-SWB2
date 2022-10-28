# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 18:38:00 2022

@author: paolo
"""

import os
import pandas as pd
os.chdir('C:/E-OBS-SWB2')

from Python.EOBSobject import EOBSobject

outpath = './Export/trials_05102022/'

inpath = './Data/E-OBS'
var = 'rr' #daily precipitation sum


# f = EOBSobject(var, inpath, outpath, folder = True, swb2 = False)
f = EOBSobject(var, inpath, outpath, folder = True, swb2 = True)
f.load()

coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)
start = 2014
end = 2018

# f.cut_spacetime(coord, start, end, saveformat='arcgrid', contourcell = 2)
f.cut_spacetime(coord, start, end, saveformat = 'netcdf', contourcell = 2)
