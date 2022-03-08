# -*- coding: utf-8 -*-
"""
Compute the total recharge

- Load the meteorological recharge generated from SWB2
- Divide the meteorological recharge in stress periods (SP). Sum the recharge
    inside the SP and save the sums in a 3D variable: rmeteo(i, nrow, ncol),
    where i is the SP number
- Load databases on land cover, irrigation recharge, municipalities of the
    area and urban recharge for each municipality
- Compute the total recharge summing the databases

@author: paolo
"""

# %% Modules

import os
import netCDF4 as nc
import numpy as np

# %% Custom functions

def leap(y):
    #input: year (int)
    #output: number of days (int)
    if((y%4 == 0) | (y%400 == 0)):
        return 366
    else:
        return 365

# %% Setup

#Paths
os.chdir('C:/E-OBS-SWB2')
swb2path = "./Model/swb2_MODELMI/output/ModelMI_net_infiltration__2014-01-01_2018-12-31__338_by_660.nc"

#Stress period definition
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 
SP4 = 107  #days, - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

# %% Load the meteorological recharge generated from SWB2

f = nc.Dataset(swb2path)

#or with a class
#f = SWB2output(swb2path)

# %% Obtain the stress periods sums

net_infiltration = np.ma.getdata(f['net_infiltration'][:,:,:])
f.close()

#outpath = "./Export/ASCII/RMeteo_SP_inch"

#Could the start and end years be retrieved from netCDF metadata?
starty = 2014 #Start year
endy = 2018   #End year

# #To save as ArcASCII, xll and yll have to be the extreme left bottom point,
# #not the center of the left-bottom cell, so remove size/2
# size = round(np.ma.getdata(f['x'][1]).item() - np.ma.getdata(f['x'][0]).item())
# xll = round(np.ma.getdata(f['x'][0]).item()) - size/2
# yll = round(np.ma.getdata(f['y'][-1]).item()) - size/2

#Extract a single year
period = range(starty, endy+1)
s, e, k = 0, 0, 0

#Creazione variabile 3D: da aggiustare, così non funziona
rmeteo = np.array(0, (len(period)*len(SPs), dim(net_infiltration)[1], dim(net_infiltration)[2]))

for y in period:
    e += leap(y)
    year = net_infiltration[s:e, :, :]
    #Extract the Stress Period
    base = 0
    for i, SP in enumerate(SPs, start = 1):
        sp = year[base:SP, :, :]
        base = SP
        #Sum the infiltration across the stress period
        #Transform into m/s
        # sp = np.sum(sp, axis = 0)*0.0254/(60*60*sp.shape[0]) #m/s
        #Keep in inches (make the transformation later). Useful because otherwise ArcGIS can't read the values (too small)
        sp = np.sum(sp, axis = 0) #inches
        
        
        #Save in the 3D variable rmeteo
        rmeteo[k, :, :] = sp
        k += 1
    s = e #It will get the subsequent day



# %% Load and organize the needed datasets



# %% Compute the total recharge


