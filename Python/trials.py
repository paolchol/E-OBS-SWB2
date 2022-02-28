# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 15:51:34 2022

@author: paolo
"""

import netCDF4 as nc
import glob
import os
import pandas as pd
import numpy as np

os.chdir('C:/E-OBS-SWB2')

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

outpath = "./Export/ASCII/Sums"
fls = glob.glob('./Data/SWB2_output/meteo/*.nc')

for i, fl in enumerate(fls, start = 1):
    f = nc.Dataset(fl)
    rainfall = np.ma.getdata(f['rainfall'][:,:,:])
    
    df = np.sum(rainfall, axis = 0)*0.0254 #meters
    df = pd.DataFrame(df)
    fname = f"{outpath}/SWB2_prec_sum_tot.asc" #To obtain a CSV just change to .csv and run this and the below line
    df.to_csv(fname, sep = ' ', header = False, index = False)
    xll = round(np.ma.getdata(f['x'][0]).item())
    yll = round(np.ma.getdata(f['y'][-1]).item())
    size = round(np.ma.getdata(f['x'][1]).item()) - xll #controlla che sia 100
    save_ArcGRID(df, fname, xll, yll, size, -9999)
    
    f.close()

#rain = np.ma.getdata(f['rainfall'][:, :, :])

# %% Sum the E-OBS precipitation over the 5 years

outpath = "./Export/ASCII/Sums"
fls = glob.glob('./Export/netCDF/calco_Daymet/*.nc')

sumtot = 0
for i, fl in enumerate(fls[0:5], start = 1):
    f = nc.Dataset(fl)
    prcp = np.ma.getdata(f['prcp'][:,:,:]) #millimeters
    
    df = np.sum(prcp, axis = 0)/1000 #meters
    xll = round(np.ma.getdata(f['x'][0]).item())
    yll = round(np.ma.getdata(f['y'][-1]).item())
    size = round(((np.ma.getdata(f['x'][1]).item() - xll + np.ma.getdata(f['y'][-2]).item() - yll))/2)
    sumtot = np.add(sumtot, df)
    f.close()

sumtot = pd.DataFrame(sumtot)
fname = f"{outpath}/EOBS_prec_sum_tot.asc" #To obtain a CSV just change to .csv and run this and the below line
sumtot.to_csv(fname, sep = ' ', header = False, index = False)
save_ArcGRID(sumtot, fname, xll, yll, size, -9999)





