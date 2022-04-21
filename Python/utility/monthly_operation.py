# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 17:18:30 2022

@author: paolo
"""

#%% Setup

import os
import pandas as pd
os.chdir('C:/E-OBS-SWB2')

from Python.EOBSobject import EOBSobject

outpath = './Export/netcdf/output_EOBSobject'

inpath = './Data/E-OBS'

coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)

start = 2014
end = 2018

# %% Run

var = ['rr', 'tn', 'tx']
#Also outpath and inpath could be provided in lists, if you want to refer to
#single files each time in different folders
for v in var:
    f = EOBSobject(inpath, v, outpath)
    f.load()
    f.cut_spacetime(coord, start, end, contourcell = 2, option = 'bundle')
    f.close_netcdf()
    
# %% Calculate monthly sum

import netCDF4 as nc
import numpy as np
import pandas as pd

def leap(y):
    #input: year (int)
    #output: number of days (int)
    if((y%4 == 0) | (y%400 == 0)):
        return 366
    else:
        return 365

def month_duration(month, year):
    if month in [1, 3, 5, 7, 8, 10, 12]:
        return 31
    elif month in [4, 6, 9, 11]:
        return 30
    else:
        if leap(year) == 365: return 28
        else: return 29

def monthly_operation(file, var, start, end, function = np.sum,
                      row = None, col = None, station = None, singlecell = False,
                      valname = 'value'):
    f = nc.Dataset(file)
    s, e, k = 0, 0, 0
    out = np.zeros((12*(end-start+1), 4))
    for year in range(start, end+1):
        for i in range(1, 12+1):
            e += month_duration(i, year)
            df = np.ma.getdata(f[var][s:e, :, :])
            op = function(df, axis = 0)
            if singlecell:
                out[k, 0] = i
                out[k, 2] = year
                out[k, 3] = op[row, col]
            k += 1
            s = e
    f.close()
    out = pd.DataFrame(out, columns = ['month', 'station', 'year', valname])
    out['station'] = station
    return out




stations = {
    'row': [4, 2, 6],
    'col': [8, 4, 10],
    'name':["Lambrate", "Busto Arsizio", "Sant'Angelo Lodigiano"]
    }
stations = pd.DataFrame(stations)

file = f"{outpath}/rr_EOBS_cut_spacetime.nc"
var = 'rr'
for _, station in stations.iterrows():
    df = monthly_operation(file, var, start, end,
                           row = station['row']-1, col = station['col']-1,
                           station = station['name'], singlecell = True,
                           valname = 'prec')
    df.to_csv(f"{outpath}/{station['name']}_prec_monthly_sum.csv", index = False)

file = f"{outpath}/tx_EOBS_cut_spacetime.nc"
var = 'tx'
for _, station in stations.iterrows():
    df = monthly_operation(file, var, start, end,
                           row = station['row']-1, col = station['col']-1,
                           station = station['name'], singlecell = True,
                           valname = 'tmax',
                           function = np.mean)
    df.to_csv(f"{outpath}/{station['name']}_tmax_monthly_mean.csv", index = False)

file = f"{outpath}/tn_EOBS_cut_spacetime.nc"
var = 'tn'
for _, station in stations.iterrows():
    df = monthly_operation(file, var, start, end,
                           row = station['row']-1, col = station['col']-1,
                           station = station['name'], singlecell = True,
                           valname = 'tmin',
                           function = np.mean)
    df.to_csv(f"{outpath}/{station['name']}_tmin_monthly_mean.csv", index = False)





