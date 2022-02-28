# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 15:05:08 2022

@author: paolo
"""

import netCDF4 as nc
import numpy as np
import os

os.chdir('C:/E-OBS-SWB2')


eobs = nc.Dataset("./Data/E-OBS/rr_ens_mean_0.1deg_reg_2011-2021_v24.0e.nc")
sample = nc.Dataset(r".\Data\tmax_Daymet_v3_2013.nc")
prova = nc.Dataset(r".\Export\netCDF\calco_Daymet\prove\prova_v7_lambert.nc")

print(sample)
print(prova)

print(sample['time'])

samplevar = sample['tmax'][0:3, :, :]
samplevar = np.ma.getdata(samplevar)

lat = sample['lat'][:]
y = sample['y'][:]

yp = prova['y'][:]
var = np.ma.getdata(prova['rr'][0:3, :, :])

print(eobs)
print(eobs['rr'])

prova.close()
sample.close()
eobs.close()


#%% Confronto tra tmax

#Samples
prec = nc.Dataset(r".\Data\prcp_Daymet_v3_2012.nc")
tmin = nc.Dataset(r".\Data\tmin_Daymet_v3_2012.nc")
tmax = nc.Dataset(r".\Data\tmax_Daymet_v3_2013.nc")
#Prove
ptmin = nc.Dataset(r".\Export\netCDF\calco_Daymet\tmin_EOBS_2014.nc")
ptmax = nc.Dataset(r".\Export\netCDF\calco_Daymet\tmax_EOBS_2014.nc")

print('prec comparison')
print(prec)
print(prec['prcp'])

print("tmin comparison")
print(tmin)
print(ptmin)
print("tmax comparison")
print(tmax)
print(tmax['tmax'])
print(ptmax['tx'])

sp = np.ma.getdata(prec['prcp'][:,:,:])
stn = np.ma.getdata(tmin['tmin'][:,:,:])
stx = np.ma.getdata(tmax['tmax'][:,:,:])
ptx = np.ma.getdata(ptmax['tmax'][:,:,:])

time = tmax['time'][:]
ptime = ptmax['time'][:]

ptmin['y'][:]

prec.close()
tmin.close()
tmax.close()
ptmin.close()
ptmax.close()

#Check the date
def x_todate(x):
    from datetime import date, timedelta
    start = date(1980,1,1)
    end = start + timedelta(days = x.item())
    #end = start + timedelta(days = x)
    return end.year, end.strftime('%m'), end.strftime('%d')

x_todate(ptmin['time'][0] - 0.5)
x_todate(tmin['time'][0] - 0.5)



