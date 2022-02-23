# -*- coding: utf-8 -*-
"""
Main script for the management of E-OBS data

@author: paolo
"""

import numpy as np
import os
import netCDF4 as nc
import glob
import pandas as pd
import utm

import time
starttime = time.time()

# %% Custom functions

def date_toeobs(y,m,d):
    #Returns how many days have passed from 1950-01-01
    from datetime import date
    start = date(1950,1,1)
    d = date(y,m,d)
    count = d - start #how many days have passed since start
    return count.days

def eobs_todate(x):
    from datetime import date, timedelta
    start = date(1950,1,1)
    end = start + timedelta(days = x.item())
    #end = start + timedelta(days = x)
    return end.year, end.strftime('%m'), end.strftime('%d')

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line + '\n' + content)
        #line.rstrip('\r\n') if you want ot remove something from line

def transfcoord(lat, lon, zoneN, zoneL):
    import utm as utm
    #The syntax is utm.from_latlon(LATITUDE, LONGITUDE)
    #The return has the form (EASTING, NORTHING, ZONE_NUMBER, ZONE_LETTER)
    #For a single value, it needs an int or a float. Change the type while giving the input
    east = []
    north = []
    if (type(lat) == float) | (type(lat) == int):
        coord = utm.from_latlon(lat, lon,
                                force_zone_number = zoneN,
                                force_zone_letter = zoneL)
        east += [coord[0]]
        north += [coord[1]]
    else:
        for i,_ in enumerate(lat):
            coord = utm.from_latlon(lat[i], lon[i],
                                force_zone_number = zoneN,
                                force_zone_letter = zoneL)
            east += [coord[0]]
            north += [coord[1]]
    return east, north

# east, north = transfcoord(lat, lon, 32, 'N')

def ArcASCIIfile(filename, ncol, nrow, xll, yll, cell, nodata):
    TheFile = open(filename,"w")
    TheFile.write("ncols 10\n")
    TheFile.write("nrows 10\n")
    TheFile.write("xllcorner     0.0\n")
    TheFile.write("yllcorner     0.0\n")
    TheFile.write("cellsize      1.0\n")
    TheFile.write("NODATA_value  -9999\n")

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

# %% Create the folders to save the exported files

os.chdir('C:/Users/paolo/Desktop/progetto E-OBS')
namefolder = ['precip','tmin','tmax']
exp = r'.\Export\ASCII'
if not os.path.exists(exp):
    os.makedirs(exp)
for name in namefolder:
    newpath = f'.\Export\ASCII\{name}'
    if not os.path.exists(newpath):
        os.makedirs(newpath)

# %% Generate daily .asc files

#Set the paths to the files
ncdf = glob.glob('C:/Users/paolo/Desktop/progetto E-OBS/Dati/E-OBS/*.nc')

#Select the needed files
# tn: minimum temperature
# tx: maximum temperature
# rr: precipitation sum

# rr = nc.Dataset(ncf[4])
# tn = nc.Dataset(ncf[6])
# tx = nc.Dataset(ncf[7])
idxs = [4, 7, 8]
tag = ['rr','tn', 'tx']
outname = ['PRCP', 'TMIN', 'TMAX']
#These tags and outnames can be put into a dataframe that can be loaded each time
#To generalize the procedure
#So that i can select the variables i want with 'rr', 'tn' ecc and than the code can
#find the idx correspondant in the file list (ncdf)

#Just to test the code
#i, idx, tt = 0, 4, 1096

for i, idx in enumerate(idxs, start = 0):
    #Load the E-OBS dataset of the variable of interest
    ds = nc.Dataset(ncdf[idx])
    #Extract longitude, latitude and time
    lon = ds['longitude'][:]
    lat = ds['latitude'][:]
    t = ds['time'][:]
    #Obtain the wanted date range using the E-OBS date format*
    start = np.where(t == date_toeobs(2014,1,1))[0][0]
    end = np.where(t == date_toeobs(2018,12,31))[0][0] + 1
    #+1 to actually consider 2018-12-31
    
    for tt in range(start, end):
        val = ds[tag[i]][tt, :, :]
        #Create a pd.DataFrame with the variable
        #Place lat and lon as column and row names
        df = pd.DataFrame(val, columns = lon, index = lat)
        df.fillna(-9999, inplace = True)
        #Cut in space using the extremes determined in the previous section
        df = df.loc[:, (df.columns > minlon) & (df.columns < maxlon)]
        df = df.loc[(df.index < maxlat) & (df.index > minlat), :]
        #Rows have to be sorted decreasing (as latitude in reality)
        #this is not true, it has been updated later!
        #df = df.sort_index(ascending = False)
        
        #Save the file
        #.asc file, with the header and format required by SWB
        #Delimiter: space
        y, m, d = eobs_todate(t[tt])
        fname = f'.\Export\ASCII\{namefolder[i]}\{outname[i]}_{y}_{m}_{d}.asc'
        #float
        round(df, 2).to_csv(fname, sep = ' ', header = False, index = False)
        #int
        #df.astype(int).to_csv(fname, sep = ' ', header = False, index = False)
        #Insert the header
        coord = utm.from_latlon(df.index[-1], df.columns[0], 32, 'N')
        #x: East, y: North
        header = f'ncols         {len(df.columns)}\nnrows         {len(df.index)}\nxllcorner     {round(coord[0])}\nyllcorner     {round(coord[1])}\ncellsize      8300\nNODATA_value  -9999'
        line_prepender(fname, header)
        
    ds.close()

#*It is the number of days passed from 1950-01-01

endtime = time.time()
print(f'Runtime of the program: {(endtime - starttime)/60}')
#1.6 min
