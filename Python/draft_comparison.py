# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 17:54:36 2022

@author: paolo
"""

# %% Load the created netCDF

f = nc.Dataset(r"C:\Users\paolo\Desktop\progetto E-OBS\Dati\Ritagli_netCDF\prova.nc", 'r')
print(f)
print(f.variables['rr'])
print(f.variables['Latitude'])
print(f.variables['Longitude'])
print(f.dimensions['time'])

f['rr'][0, :, :]
f.close()

# %% Check the sample dataset provided by SWB2 manual

sample = nc.Dataset(r"C:\Users\paolo\Desktop\progetto E-OBS\Dati\prcp_Daymet_v3_2012.nc")

print(sample)
#Dimensions
print(sample.dimensions['x'])
print(sample.dimensions['time'])

#Variables
print(sample.variables['time'])
print(sample['time_bnds'])
print(sample['prcp'])
print(sample['lambert_conformal_conic'])
print(sample.variables['x'])
print(sample['yearday'])
print(sample['lat'])
print(sample['lon'])

time = sample['time'][:]
time_bnds = sample['time_bnds'][:]
x = sample['x'][:]
y = sample['y'][:]
lat = sample['lat'][:]
lon = sample['lon'][:]
yearday = sample['yearday'][:]

lcc = sample['lambert_conformal_conic'][:]
lcc = np.ma.getdata(lcc)

sample.close()


print(ncf)


def date_toeobs(y,m,d):
    #Returns how many days have passed from 1950-01-01
    from datetime import date
    start = date(1950,1,1)
    d = date(y,m,d)
    count = d - start #how many days have passed since start
    return count.days

date_toeobs(2012, 1, 1)

def eobs_todaymet(y):
    from datetime import date
    dstart = date(1980, 1, 1)
    estart = date(1950, 1, 1)
    k = dstart - estart
    y1 = y - k.days + 0.5
    return y1

eobs_todaymet(date_toeobs(2012, 1, 1))

eobs_todaymet(yy)


