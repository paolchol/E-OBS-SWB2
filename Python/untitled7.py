# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 15:43:17 2022

@author: AMARANTO
"""

import netCDF4 as nc
import glob
import os
import pandas as pd
import numpy as np

os.chdir('C:/Users/paolo/Desktop/progetto E-OBS/Dati/E-OBS')
fls = glob.glob('*.nc')

ncf = nc.Dataset(fls[6]) #tg
la = ncf['latitude'][:]
lo = ncf['longitude'][:]
# crei idx per lat e lon

"""
Esempio
lmin = 33
lmax = 48
idx_lat = np.intersect1d(np.where(la > lmin)    , np.where(la < lmax)  )


"""


dt = pd.date_range(start='1950-01-01', end='2021-06-30')
yR = dt.year
yU = yR.unique()

y = 2014

for y in yU:
    
    idx = np.where(yR == y)
    wi = ncf['tg'][idx[0], :, :] # metti qua idx la e lon
    
    fl = str(y) + '_tg_orig.nc'
    ds = nc.Dataset(fl, 'w', format = 'NETCDF4')
    
    lat = ds.createDimension('lat', len(la)) # aggiorna dimensioni spaziali
    lon = ds.createDimension('lon', len(lo))
    tim = ds.createDimension('tim', wi.shape[0])

    #times = ds.createVariable('timed', 'f4', ('timed',))
    latitude = ds.createVariable('lat', 'f4', ('lat',))
    longitude = ds.createVariable('lon', 'f4', ('lon',))
    time = ds.createVariable('tim', 'f4', ('tim',))
    value = ds.createVariable('tg', 'f4', ('tim','lat', 'lon',))
    value.units = 'nonmiricordo'
    
    latitude[:] = la
    longitude[:] = lo
    time[:] =  dt[idx[0]]
    
    value[:] = wi
    
    ds.close()

ncf.close()    
