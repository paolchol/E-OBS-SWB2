# -*- coding: utf-8 -*-

import os
import pandas as pd
os.chdir('C:/E-OBS-SWB2')

from Python.EOBSobject import EOBSobject

outpath = './Export'
inpath = './Data/E-OBS'
var = 'rr'
f = EOBSobject(var, inpath, outpath, folder = True)
f.load()

SPs = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)
start = 2014
end = 2018

f.SP_sum(SPs,    coord = coord, start = start, end = end,
                 contourcell = 2,
                 method = 'cut_spacetime',
                 units = 'mm')
f.close_netcdf()