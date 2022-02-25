# -*- coding: utf-8 -*-
"""
Creation of netCDF files compatible with SWB2 software
- Starting point: E-OBS data
- netCDF building based on the Daymet netCDF provided by SWB2 manual

@author: paolo
"""

# %% Modules

import netCDF4 as nc
import glob
import os
import pandas as pd
import numpy as np

# %% Project's directory
os.chdir('C:/E-OBS-SWB2')

# %% Custom functions

def transf(lat_t, lon_t, zoneN, zoneL, var = 'xy'):
    #Transforms the lat/lon data provided in x/y projected UTM coordinates
    import utm
    x = []
    y = []
    for i in lat_t:
        for j in lon_t:
            #The syntax is utm.from_latlon(LATITUDE, LONGITUDE)
            #The return has the form (EASTING, NORTHING, ZONE_NUMBER, ZONE_LETTER)
            xj, yi, _, _ = utm.from_latlon(i, j,
                            force_zone_number = zoneN,
                            force_zone_letter = zoneL)
            x += [xj]
        y += [yi]
    x = x[0:len(lon_t)]
    if(var == 'xy'):
        return x, y
    elif(var == 'x'):
        return x
    elif(var == 'y'):
        return y

def eobs_todaymet(y):
    #redefine: starting from 1980, adding 0.5
    from datetime import date
    dstart = date(1980, 1, 1)
    estart = date(1950, 1, 1)
    k = dstart - estart
    y1 = y - k.days + 0.5
    return y1

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

# %% Set the parameters

index = [4, 7, 8] # referred to the index of the files in fls = glob.glob('./Dati/E-OBS/*.nc')
tag = ['rr','tn', 'tx'] #names in the E-OBS files
outname = ['prcp', 'tmin', 'tmax'] #names in the output files (Daymet copies)
units = ['mm/day', 'degrees C', 'degrees C']

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

#np.where(t == '2016-02-29')
# %% Generate the new files

#Get the path to the files
outpath = r'./Export/netCDF/calco_Daymet'
#outpathswb2 = r'./swb2_MODELMI/climate_ncfile'
fls = glob.glob('./Data/E-OBS/*.nc')

# n = 4
# year = 2014
# i = 2

for i, n in enumerate(index, start = 0):
    ncf = nc.Dataset(fls[n])
    la = ncf['latitude'][:]
    lo = ncf['longitude'][:]
    
    idx_lat = np.intersect1d(np.where(la > minlat), np.where(la < maxlat))
    idx_lon = np.intersect1d(np.where(lo > minlon), np.where(lo < maxlon))
        
    lox = transf(la[idx_lat], lo[idx_lon], 32, 'N', 'x')
    lay = transf(la[idx_lat], lo[idx_lon], 32, 'N', 'y')
    
    #Sort latitude descending
    Mlat = pd.DataFrame(la[idx_lat]).sort_values(0, ascending = False)
    Mlat = pd.concat([Mlat]*len(idx_lon), axis = 1, ignore_index = True)
    Mlat.index = range(0,len(idx_lat))
    Mlon = pd.DataFrame(columns = list(range(0,len(idx_lon))))
    Mlon.loc[0] = lo[idx_lon]
    Mlon = pd.concat([Mlon]*len(idx_lat), axis = 0, ignore_index = True)
    
    for year in yU[np.where((yU >= 2014) & (yU <= 2018))]:
        
        #idx = np.where((yR == year) & (dt != '2016-02-29'))
        idx = np.where(yR == year)
        val = ncf[tag[i]][idx[0], idx_lat, idx_lon]
        val = np.ma.getdata(val)
        
        #To be coherent with the convention used in SWB for data sorting
        #Flip the variable matrix around the horizontal axis
        val = np.flip(val, axis = 1)
        #Sort the y coordinate in descending order
        lay.sort(reverse = True)
        
        fname = f'{outpath}/{outname[i]}_EOBS_{year}.nc'
        #fname = f'{outpath}/prove/prova_v7_lambert.nc'
        ds = nc.Dataset(fname, 'w', format = "NETCDF3_CLASSIC")
        
        ##General metadata
        ds.description = f"E-OBS dataset over the Ticino-Adda groundwater basin, year {year}"
        ds.source = "E-OBS v24.0"
        ds.start_day = f"01/01/{year}"
        ds.author = "paolocolombo1996@gmail.com"
        
        ## Dimensions
        ds.createDimension('x', len(idx_lon))
        ds.createDimension('y', len(idx_lat))
        ds.createDimension('time', None)
        ds.createDimension('nv', 2)
        
        ## Variables
        #lambert_conformal_conic()
        lcc = ds.createVariable('lambert_conformal_conic', 'h')
        lcc.description = 'fake variable just to try a point, this data was not obtained with a lambert conformal conic'
        #x(x)
        x = ds.createVariable('x', 'd', ('x'))
        x.units = 'm'
        x.long_name = 'x coordinate of projection'
        x.standard_name = 'projection_x_coordinate'
        #y(y)
        y = ds.createVariable('y', 'd', ('y'))
        y.units = 'm'
        y.long_name = 'y coordinate of projection'
        y.standard_name = 'projection_y_coordinate'
        #Time
        time = ds.createVariable('time', 'd', ('time'))
        time.units = 'days since 1980-01-01 00:00:00 UTC'
        time.calendar = '2016 leap day removed'
        time.bounds = 'time_bnds'
        #Latitude
        lat = ds.createVariable('lat', 'd', ('y','x'))
        lat.units = 'degrees_north'
        lat.long_name = 'latitutde coordinate'
        lat.standard_name = 'latitude'
        #Longitude
        lon = ds.createVariable('lon', 'd', ('y','x'))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude coordinate'
        lon.standard_name = 'longitude'
        #Day of the year
        yearday = ds.createVariable('yearday', 'h', ('time'))
        #Value
        value = ds.createVariable(outname[i], 'f4', ('time','y','x'), fill_value = -9999)
        value.units = units[i]
        value.missing_value = -9999.0
        value.coordinates = 'lat lon'
        #Time bounds
        time_bnds = ds.createVariable('time_bnds', 'd', ('time','nv'))
        
        ## Fill the variables
        x[:] = lox
        y[:] = lay
        lat[:] = Mlat
        lon[:] = Mlon
        tout = eobs_todaymet(np.where(tyR == year)[0])
        time[:] = tout
        yearday[:] = range(0,365)
        time_bnds[:] = pd.DataFrame({'0':tout-0.5, '1':tout+0.5})
        value[:] = np.around(val, 1)
        
        #Close the file
        ds.close()
    print(f'Variable: {tag[i]} (E-OBS), {outname[i]} (Daymet)')
    print(f'Number of rows: {len(idx_lat)}\nNumber of columns: {len(idx_lon)}')
    print(f'Latitude (Y) and Longitude (X) of extention:\nX0: {lo[idx_lon][0]}\nY0: {la[idx_lat][0]}\nX1: {lo[idx_lon][-1]}\nY1: {la[idx_lat][-1]}')
    
    ncf.close()
