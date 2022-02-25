# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 17:48:39 2022

@author: paolo
"""

import netCDF4 as nc
import glob
import os
import pandas as pd
import numpy as np

# %% Define area of the project

#Points obtained from the project's perimeter*
points = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
model_extremes = pd.DataFrame(points)
#*this can be updated to read directly the point shapefile

#Take min and max of the coordinates to define the angles of the desired rectangle
#0.1 (degrees) is added to the maxs and subtracted to the mins to take a cell
#more to interpolate later without border effects
minlon = min(model_extremes.lon) - 0.1
maxlon = max(model_extremes.lon) + 0.1
minlat = min(model_extremes.lat) - 0.1
maxlat = max(model_extremes.lat) + 0.1

# %% Set the parameters

index = [4, 7, 8]
tag = ['rr','tn', 'tx']
outname = ['prcp', 'tmin', 'tmax']
units = ['mm', 'C degrees', 'C degrees']

# %% Get time indexes

#Insert start and end dates of E-OBS
#In the files we have, the E-OBS datasets were cut between 2011 and 2021
#The date inside the variable 'time' is still in the original format,
#registered as the number of days passed from 1950-01-01
#The indexes dt and yR however will be used to find the position in the array, not
#the actual relative date
dt = pd.date_range(start = '2011-01-01', end = '2021-06-30')
yR = dt.year
yU = yR.unique()

#E-OBS date definition (to insert it in the output .nc file)
t = pd.date_range(start = '1950-01-01', end = '2021-06-30')
tyR = t.year

# %% Generate the new files

#Get the path to the files
os.chdir('C:/Users/paolo/Desktop/progetto E-OBS/Dati/E-OBS')
outpath = 'C:/Users/paolo/Desktop/progetto E-OBS/Dati/Ritagli_netCDF'
fls = glob.glob('*.nc')

n = 4
y = 2014
i = 0

for i, n in enumerate(index, start = 0):
    ncf = nc.Dataset(fls[n])
    la = ncf['latitude'][:]
    lo = ncf['longitude'][:]
    
    idx_lat = np.intersect1d(np.where(la > minlat), np.where(la < maxlat))
    idx_lon = np.intersect1d(np.where(lo > minlon), np.where(lo < maxlon))
    
    for y in yU[np.where((yU >= 2014) & (yU <= 2018))]:
        
        idx = np.where(yR == y)
        val = ncf[tag[i]][idx[0], idx_lat, idx_lon]
        val = np.ma.getdata(val)
        
        #fname = f'{outpath}/{outname[i]}_EOBS_{str(y)}.nc'
        fname = f'{outpath}/prova_v2.nc'
        ds = nc.Dataset(fname, 'w', format = 'NETCDF4')
        
        ds.createDimension('lat', len(idx_lat))
        ds.createDimension('lon', len(idx_lon))
        ds.createDimension('time', val.shape[0])
        
        latitude = ds.createVariable('Latitude', 'd', ('lat'))
        longitude = ds.createVariable('Longitude', 'f4', ('lon'))
        value = ds.createVariable(tag[i], 'f4', ('time','lat', 'lon'), fill_value = -9999)
        time = ds.createVariable('Time', 'i4', ('time'))
        time.units = 'days since 1950-01-01 00:00:00 UTC'
        value.units = units[i]
        
        tout = np.where(tyR == y)
        latitude[:] = np.float32(np.ma.getdata(la[idx_lat]))
        longitude[:] = np.float32(np.ma.getdata(lo[idx_lon]))
        time[:] = np.int32(tout[0])
        value[:] = val
        
        ds.close()
    
    ncf.close()
