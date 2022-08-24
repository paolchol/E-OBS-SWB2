# -*- coding: utf-8 -*-
"""
Creation of ArcASCII GRID files

- Starting point: E-OBS data
- ArcASCII GRID files created following ArcASCII standards

The goal was to create these files to use them as input for SWB1 and SWB2.
The generated ArcASCII are anyway in the standard format, and thus applicable to generate
this type of file for any use.

@author: paolo
"""

# %% Modules

import numpy as np
import os
import netCDF4 as nc
import glob
import pandas as pd
import utm
import time
import math

#%% Project directory
os.chdir('C:/E-OBS-SWB2')

# %% Custom functions

from Python.custom_functions import date_toeobs
from Python.custom_functions import eobs_todate
from Python.custom_functions import save_ArcGRID

# %% Define area of the project

#Points obtained from the project's perimeter*
points = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
model_extremes = pd.DataFrame(points)
#*this can be updated to read directly the point shapefile
#Or even directly the area shapefile, then search for the maximum and minimum points

#Take min and max of the coordinates to define the angles of the desired rectangle
#0.1 (degrees) is added to the maxs and subtracted to the mins to take a cell
#more to interpolate later without border effects
#an additional 0.1 (degrees) are used to cover the "ancillary" points in SWB
minlon = min(model_extremes.lon) - 0.1 - 0.1
minlat = min(model_extremes.lat) - 0.1 - 0.1
maxlon = max(model_extremes.lon) + 0.1 + 0.1
maxlat = max(model_extremes.lat) + 0.1 + 0.1

#Ancillary points of SWB directly loaded and considered as starting point to cut
# f = pd.read_csv("./Data/ancillary_grid_extentions_SWB.txt", index_col = 0, sep = "\t")
# minlat = f['y'][0] - 0.1
# minlon = f['x'][0] - 0.1
# maxlat = f['y'][1] + 0.1
# maxlon = f['x'][1] + 0.1

# %% Create the folders to save the exported files

#os.chdir('C:/E-OBS-SWB2')
namefolder = ['precip','tmin','tmax']
exp = r'.\Export\ASCII'
if not os.path.exists(exp):
    os.makedirs(exp)
for name in namefolder:
    newpath = f'.\Export\ASCII\{name}'
    if not os.path.exists(newpath):
        os.makedirs(newpath)

# %% Set the parameters

#Select the needed files
# tn: minimum temperature
# tx: maximum temperature
# rr: precipitation sum

# rr = nc.Dataset(fls[4])
# tn = nc.Dataset(fls[6])
# tx = nc.Dataset(fls[8])
idxs = [4, 7, 8]
tag = ['rr','tn', 'tx']
outname = ['PRCP', 'TMIN', 'TMAX']
#These tags and outnames can be put into a dataframe that can be loaded each time
#To generalize the procedure
#So that i can select the variables i want with 'rr', 'tn' ecc and than the code can
#find the idx correspondant in the file list (fls)

#Just to test the code
#i, idx, tt, j = 0, 4, 1096, 1

# %% Generate daily .asc files

starttime = time.time()

#Set the paths to the files
fls = glob.glob('./Data/E-OBS/*.nc')

for i, idx in enumerate(idxs, start = 0):
    #Load the E-OBS dataset of the variable of interest
    ncf = nc.Dataset(fls[idx])
    #Extract longitude, latitude and time from E-OBS
    lon = ncf['longitude'][:]
    lat = ncf['latitude'][:]
    t = ncf['time'][:]
    #Obtain the wanted date range using the E-OBS date format
    start = np.where(t == date_toeobs(2014,1,1))[0][0]
    end = np.where(t == date_toeobs(2018,12,31))[0][0] + 1 #+1 to actually consider 2018-12-31
    #Obtain the range of lat and lon to extract from the whole dataset
    idx_lat = np.intersect1d(np.where(lat > minlat), np.where(lat < maxlat))
    idx_lon = np.intersect1d(np.where(lon > minlon), np.where(lon < maxlon))
    #Extract the values in the time and space range
    val = ncf[tag[i]][range(start, end), idx_lat, idx_lon]
    #Flip the variable matrix around the horizontal axis
    val = np.flip(val, axis = 1)
    
    size = round(lat[idx_lat][1] - lat[idx_lat][0], 1)
    xll = lon[idx_lon][0] - size/2
    yll = lat[idx_lat][0] - size/2
    
    # xll, yll, _, _ = utm.from_latlon(lat[idx_lat][0], lon[idx_lon][0], 32, 'N')
    # x, y, _, _ = utm.from_latlon(lat[idx_lat][1], lon[idx_lon][1], 32, 'N')
    # size = round(math.sqrt((x-xll)*(y-yll))) #side of a square with the same area as the rectangle
    
    #Cycle through the days
    for j, tt in enumerate(range(start, end)):
        #Extract the corresponding day
        df = pd.DataFrame(val[j, :, :])
        
        #Save the dataset as an .asc file
        #The header and format required by SWB (defined as Arc ASCII GRID)
        #Delimiter: space
        y, m, d = eobs_todate(t[tt], number = True)
        fname = f'.\Export\ASCII\{namefolder[i]}\{outname[i]}_{y}_{m}_{d}.asc'
        save_ArcGRID(round(df, 1), fname, xll, yll, size, -9999)
        
    ncf.close()

endtime = time.time()
print(f'Runtime of the program: {(endtime - starttime)/60}')
#1.14 min
