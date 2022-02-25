# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 14:43:39 2022

@author: paolo
"""

#Make up the data
import netCDF4 as nc
import os
import numpy as np

os.chdir('C:/Users/paolo/Desktop/progetto E-OBS/Dati/E-OBS/prove')

lon = np.arange(45,101,2)
lat = np.arange(-30,25,2.5)
z = np.arange(0,200,10)
x = np.random.randint(10,25, size=(len(lon), len(lat), len(z)))
noise = np.random.rand(len(lon), len(lat), len(z))
temp_data = x+noise

#Generate the dataset
f = nc.Dataset('sample_v2.nc', 'w', format = 'NETCDF4') #'w' stands for write
#A group can be created: directory/folder within the dataset
#tempgrp = f.createGroup('Temp_data')
#With a group, the next steps (until #Add global attributes) would be all doene with tempgroup.function

#Specify the dimensions
f.createDimension('lon', len(lon))
f.createDimension('lat', len(lat))
f.createDimension('z', len(z))
f.createDimension('time', None)

#Building variables
longitude = f.createVariable('Longitude', 'f4', 'lon')
latitude = f.createVariable('Latitude', 'f4', 'lat')
levels = f.createVariable('Levels', 'i4', 'z')
temp = f.createVariable('Temperature', 'f4', ('time', 'lat', 'lon', 'z'))
time = f.createVariable('Time', 'i4', 'time')

#Passing data into variables
longitude[:] = lon #The "[:]" at the end of the variable instance is necessary
latitude[:] = lat
levels[:] = z
temp[0,:,:,:] = temp_data

#get time in days since Jan 01,01
from datetime import datetime
today = datetime.today()
time_num = today.toordinal()
time[0] = time_num

#Add global attributes
f.description = "Example dataset containing one group"
f.history = "Created " + today.strftime("%d/%m/%y")

#Add local attributes to variable instanceslongitude.units = 'degrees east'
latitude.units = 'degrees north'
time.units = 'days since Jan 01, 0001'
temp.units = 'Kelvin'
levels.units = 'meters'
temp.warning = 'This data is not real!'

#Look what has been done
print(f)

#Final step
f.close()