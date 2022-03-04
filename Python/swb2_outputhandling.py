# -*- coding: utf-8 -*-
"""
SWB2 output handling

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

def leap(y):
    if((y%4 == 0) | (y%400 == 0)):
        d = 366
    else:
        d = 365
    return d

def save_ArcGRID(df, fname, xll, yll, size, nodata):
    #df has to be a Pandas DataFrame
    #xll: x coordinate of the left bottom corner (lon)
    #yll: y coordinate of the left bottom corner (lat)
    #size: cell size (m)
    #nodata: value assigned to nodata
    def line_prepender(filename, line):
        with open(filename, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line + '\n' + content)
            #line.rstrip('\r\n') if you want ot remove something from line
    df.to_csv(fname, sep = ' ', header = False, index = False)
    header = f'ncols         {len(df.columns)}\nnrows         {len(df.index)}\nxllcorner     {xll}\nyllcorner     {yll}\ncellsize      {size}\nNODATA_value  {nodata}'
    line_prepender(fname, header)

# %% Create the sum over 5 years for different tests

outpath = "./Export/ASCII/RMeteo_tot" #not include /
inpath = "./Data/SWB2_output/" #include /
fls = glob.glob(f'{inpath}*.nc')

variable = 'net_infiltration'

for i, fl in enumerate(fls, start = 1):    
    f = nc.Dataset(fl)
    #net_infiltration = np.ma.getdata(f['net_infiltration'][:,:,:])
    
    df = np.sum(np.ma.getdata(f[variable][:,:,:]), axis = 0)*0.0254 #meters
    df = pd.DataFrame(df)
    
    size = round(np.ma.getdata(f['x'][1]).item() - np.ma.getdata(f['x'][0]).item()) #controlla che sia 100
    xll = round(np.ma.getdata(f['x'][0]).item()) - size/2
    yll = round(np.ma.getdata(f['y'][-1]).item()) - size/2
    
    fname = f"{outpath}/{fl[len(inpath):].split('_', 1)[0]}_{variable}_sum_tot.asc" #To obtain a CSV just change to .csv and run this and the below line
    save_ArcGRID(df, fname, xll, yll, size, -9999)
    
    f.close()

# %% Get the sum over the 4 stress periods

# swbout = "./swb2_MODELMI/output"
# f = nc.Dataset(f'{swbout}/ModelMI_net_infiltration__2014-01-01_2018-12-31__338_by_660.nc')

fls = glob.glob('./Data/SWB2_output/*.nc')
f = nc.Dataset(fls[0]) #here update searching for "net_infiltration" in the filename
net_infiltration = np.ma.getdata(f['net_infiltration'][:,:,:])

outpath = "./Export/ASCII/RMeteo_SP_inch"

#Stress period definition
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 
SP4 = 107  #days, - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

#print(f['time'])
#units: days since 2014-01-01
#it includes the leap year in 2016! It has to be considered

#Could the start and end years be retrieved from netCDF metadata?
starty = 2014 #Start year
endy = 2018   #End year

#To save as ArcASCII, xll and yll have to be the extreme left bottom point,
#not the center of the left-bottom cell, so remove size/2
size = round(np.ma.getdata(f['x'][1]).item() - np.ma.getdata(f['x'][0]).item())
xll = round(np.ma.getdata(f['x'][0]).item()) - size/2
yll = round(np.ma.getdata(f['y'][-1]).item()) - size/2

#Extract a single year
period = range(starty, endy+1)
s = 0
e = 0
for y in period:
    e += leap(y)
    year = net_infiltration[s:e, :, :]
    #Extract the Stress Period
    base = 0
    for i, SP in enumerate(SPs, start = 1):
        #print(i, SP)
        sp = year[base:SP, :, :]
        base = SP
        #Sum the infiltration across the stress period
        #Transform into m/s
        #sp = np.sum(sp, axis = 0)*0.0254/(60*60*sp.shape[0])
        #Keep in inches (make the transformation later). Useful because otherwise ArcGIS can't read the values (too small)
        sp = np.sum(sp, axis = 0)
        #Save as Arc GRID ASCII file
        fname = f'{outpath}/RMeteo_{y}_SP{i}_inch.asc'
        sp = pd.DataFrame(sp)
        save_ArcGRID(sp, fname, xll, yll, size, -9.9e-20)
        #-9.9e-20 got from _FillValue of 'net_infiltration' in the netCDF file
        #print(sp.shape)
    #print(s, e, year.shape)
    s = e #It will get the subsequent day

f.close()
