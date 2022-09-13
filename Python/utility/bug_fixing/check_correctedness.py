# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 17:32:31 2022

@author: paolo
"""

import os
import pandas as pd
import numpy as np
import netCDF4 as nc

from Python.EOBSobject import EOBSobject

#%% Generate the trials
os.chdir('C:/E-OBS-SWB2')
outpath = './Export/netCDF'
inpath = './Data/E-OBS'
var = 'rr' #daily precipitation sum
outname = 'prova1309'
start = 2018
end = 2019

#standard
f = EOBSobject(var, inpath, outpath, outname = f'{outname}', folder = True)
f.load()

#cut in space
coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)
f.cut_space(coord)
f.cut_spacetime(coord, start, end)
f.cut_time(2017, 2017)

#SWB2
f = EOBSobject(var, inpath, outpath, outname = f'{outname}_swb2', folder = True,
               swb2=True)
f.load()

#cut in space and time
coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)
f.cut_space(coord)
f.cut_spacetime(coord, start, end)
f.cut_time(2017, 2017)

#%% Make comparisons

cut = nc.Dataset(f"{outpath}/prova1309_EOBS_cut_space.nc")
cut_10 = np.ma.getdata(cut['prova1309'][1:10, :, :])
res = f.cut_space(coord, autosave = False, ext = True)
o_10 = np.ma.getdata(f.netcdf['rr'][1:10, res['idx_lat'], res['idx_lon']])


#confrontare i vettori di latitudine e longitudine dei due netcdf

# %% Get the dates in e-obs and daymet formats

from Python.custom_functions import date_toeobs
date_toeobs(2017, 1, 31)

from Python.custom_functions import date_fromstart
from datetime import date
date_fromstart(2017, 1, 31, start = date(1980, 1, 1))
