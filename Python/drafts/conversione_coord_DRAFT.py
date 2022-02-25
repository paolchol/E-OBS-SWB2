# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 12:58:25 2022

@author: paolo
"""



lat_t = np.ma.getdata(la[idx_lat])
lon_t = np.ma.getdata(lo[idx_lon])

lat_t = la[idx_lat]
lon_t = lo[idx_lon]

zoneN = 32
zoneL = 'N'


def transf(lat_t, lon_t, zoneN, zoneL, var = 'xy'):
    #The syntax is utm.from_latlon(LATITUDE, LONGITUDE)
    #The return has the form (EASTING, NORTHING, ZONE_NUMBER, ZONE_LETTER)
    from utm import utm
    x = []
    y = []
    for i in lat_t:
        for j in lon_t:
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

x = transf(lat_t, lon_t, 32, 'N', 'x')
