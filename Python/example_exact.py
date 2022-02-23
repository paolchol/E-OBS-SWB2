# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 17:10:52 2022

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
f = nc.Dataset('sample_exact.nc', 'w', format = 'NETCDF4') #'w' stands for write
#A group can be created: directory/folder within the dataset
tempgrp = f.createGroup('Temp_data')

#Specify the dimensions
tempgrp.createDimension('lon', len(lon))
tempgrp.createDimension('lat', len(lat))
tempgrp.createDimension('z', len(z))
tempgrp.createDimension('time', None)

#Building variables
longitude = tempgrp.createVariable('Longitude', 'f4', 'lon')
latitude = tempgrp.createVariable('Latitude', 'f4', 'lat')  
levels = tempgrp.createVariable('Levels', 'i4', 'z')
temp = tempgrp.createVariable('Temperature', 'f4', ('time', 'lon', 'lat', 'z'))
time = tempgrp.createVariable('Time', 'i4', 'time')

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

#Add local attributes to variable instances
longitude.units = 'degrees east'
latitude.units = 'degrees north'
time.units = 'days since Jan 01, 0001'
temp.units = 'Kelvin'
levels.units = 'meters'
temp.warning = 'This data is not real!'

#Look what has been done
print(f)

#Final step
f.close()