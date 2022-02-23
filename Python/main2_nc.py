# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 11:40:00 2022

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

index = [5, ]
tag = ['rr','tn', 'tx']
outname = ['prcp', 'tmin', 'tmax']

# %% Operations on netCDF

os.chdir('C:/Users/paolo/Desktop/progetto E-OBS/Dati/E-OBS')
fls = glob.glob('*.nc')

ncf = nc.Dataset(fls[5]) #rr
#print(ncf)

la = ncf['latitude'][:]
lo = ncf['longitude'][:]

idx_lat = np.intersect1d(np.where(la > minlat), np.where(la < maxlat))
idx_lon = np.intersect1d(np.where(lo > minlon), np.where(lo < maxlon))

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

os.chdir('C:/Users/paolo/Desktop/progetto E-OBS/Dati/Ritagli_netCDF')

y = 2014

#for y in yU[np.where((yU >= 2014) & (yU <= 2018))]:
    
idx = np.where(yR == y)
val = ncf.variables['rr'][idx[0], idx_lat, idx_lon]

fname = f'prcp_EOBS_{str(y)}_v7.nc'
ds = nc.Dataset(fname, 'w', format = 'NETCDF4')

ds.createDimension('lat', len(idx_lat))
ds.createDimension('lon', len(idx_lon))
ds.createDimension('time', val.shape[0])

latitude = ds.createVariable('Latitude', 'f4', ('lat'))
longitude = ds.createVariable('Longitude', 'f4', ('lon'))
value = ds.createVariable('rr', 'f4', ('time','lat', 'lon'))
time = ds.createVariable('Time', 'i4', ('time'))
value.units = 'inserireunit'

tout = np.where(tyR == y)
latitude[:] = np.float32(np.ma.getdata(la[idx_lat]))
longitude[:] = np.float32(np.ma.getdata(lo[idx_lon]))
time[:] = np.int32(tout[0])
value[:] = np.ma.getdata(val)
#Check how the arrays are given to the variables

#Fill 


#check the date format needed for SWB

ds.close()

ncf.close()    

# %% Load the created netCDF

f = nc.Dataset(fname, 'r')
print(f)
print(f.variables['rr'])
print(f.variables['Latitude'])
print(f.variables['Longitude'])

f['rr'][0, :, :]
