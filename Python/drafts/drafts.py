# -*- coding: utf-8 -*-
"""
Drafts
Created on Mon Jan 31 16:59:38 2022

@author: paolo
"""

# %% Drafts


#32632

# import shapefile
# points = shapefile.Reader("./Dati/estremi_modello_wgs84.shp")

# pp = points.records()
# pp[3].lat

# dic = pp[3].as_dict()

# feature = points.shapeRecords()
# first = feature.shape.__geo_interface__

#ora copia incollo, poi lo automatizzerò

model_extremes = [(8.691, 45.611),
                  (8.929, 45.308),
                  (9.524, 45.610),
                  (9.537, 45.306)]

max(model_extremes)


#interpolate, create a 100x100 gridded file: create a separate function
#put a condition on the tag to change the unit measure

lon = ds['longitude'][:]
lat = ds['latitude'][:]
t = ds['time'][:]

t_mean = ds['tg'][:]

#dt = pd.date_range(start="1950-01-01",end="2021-06-30")

#idx1 = np.where()

dt = pd.date_range(start="2014-01-01",end="2018-12-31")


for i in range(0, len(t)):
    
    t_mean = ds['tg'][i, :, :]
    fnam = 'temperature_' + str(i) + '_.txt'
    np.savetxt(fnam, t_mean)

ds.close()

with open(fname, 'a') as f:
    f.write(header + "\n")
    f.close