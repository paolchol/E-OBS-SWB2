# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 17:38:10 2022

@author: AMARANTO
"""

import numpy as np
import sys
import os
import netCDF4 as nc
import glob
import time
import pandas as pd
import rpy2


ncf = glob.glob('*.nc')
ds = nc.Dataset(ncf[0])

lon = ds['longitude'][:]
lat = ds['latitude'][:]
t = ds['time'][:]

t_mean = ds['tg'][:]

#dt = pd.date_range(start="1950-01-01",end="2021-06-30")

#idx1 = np.where()

dt = pd.date_range(start="2014-01-01",end="2018-12-31")


for i in range(0, len(t)):
    
    t_mean = ds['tg'][i, :, :]
    
    
    fnam = 'temperature_' + str(i) + '_.txt'
    np.savetxt(fnam, t_mean)



    
    
ds.close()