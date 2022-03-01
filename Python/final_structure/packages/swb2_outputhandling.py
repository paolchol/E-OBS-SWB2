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

# %% Open and get net infiltration

# swbout = "./swb2_MODELMI/output"
# f = nc.Dataset(f'{swbout}/ModelMI_net_infiltration__2014-01-01_2018-12-31__338_by_660.nc')

# print(f)
# print(f['net_infiltration'])

# net_infiltration = np.ma.getdata(f['net_infiltration'][:,:,:])

# Create the sum for different tests
outpath = "./Export/ASCII/RMeteo_tot"
fls = glob.glob('./Data/SWB2_output/*.nc')

for i, fl in enumerate(fls, start = 1):
    f = nc.Dataset(fl)
    net_infiltration = np.ma.getdata(f['net_infiltration'][:,:,:])
    
    df = np.sum(net_infiltration, axis = 0)*0.0254 #meters
    df = pd.DataFrame(df)
    fname = f"{outpath}/p{i}_sum_tot.asc" #To obtain a CSV just change to .csv and run this and the below line
    df.to_csv(fname, sep = ' ', header = False, index = False)
    xll = round(np.ma.getdata(f['x'][0]).item())
    yll = round(np.ma.getdata(f['y'][-1]).item())
    size = round(np.ma.getdata(f['x'][1]).item()) - xll #controlla che sia 100
    save_ArcGRID(df, fname, xll, yll, size, -9999)
    
    f.close()

# %% Get a sum over the 5 years to compare with the results of SWB1

df = np.sum(net_infiltration, axis = 0)*0.0254 #meters
df = pd.DataFrame(df)
fname = "somma_5anni.asc" #To obtain a CSV just change to .csv and run this and the below line
df.to_csv(fname, sep = ' ', header = False, index = False)

xll = round(np.ma.getdata(f['x'][0]).item())
yll = round(np.ma.getdata(f['y'][-1]).item())
size = size = round(np.ma.getdata(f['x'][1]).item()) - xll
# header = f'ncols         {len(df.columns)}\nnrows         {len(df.index)}\nxllcorner     {xll}\nyllcorner     {yll}\ncellsize      {size}\nNODATA_value  -9999'
# line_prepender(fname, header)
save_ArcGRID(df, fname, xll, yll, size, -9999)

# %% Get the sum over the 4 stress periods

outpath = "./Export/ASCII/RMeteo_SP_inch"

#Stress period definition
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 
SP4 = 107  #days, - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

print(f['time'])
#units: days since 2014-01-01
#it includes the leap year in 2016! It has to be considered

#Could be retrieved from netCDF metadata?
starty = 2014 #Start year
endy = 2018   #End year

xll = round(np.ma.getdata(f['x'][0]).item())
yll = round(np.ma.getdata(f['y'][-1]).item())
size = round(np.ma.getdata(f['x'][1]).item()) - xll
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

#%% Trials
def eobs_todate(x):
    from datetime import date, timedelta
    start = date(2014,1,1)
    #end = start + timedelta(days = x.item())
    end = start + timedelta(days = x)
    return end.year, end.strftime('%m'), end.strftime('%d')

eobs_todate(1461)



