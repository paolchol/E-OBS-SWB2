# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 16:27:24 2022

@author: paolo
"""


import pandas as pd
import numpy as np

from Python.RechargeCalc import RechargeCalc

startyear = 2014
endyear = 2018
cell_area = 100*100 #m2
swb2path = "./Data/SWB2_output/VersioneFINALE_net_infiltration.nc"
inputpath = "./Data/Calcolo_ricarica_totale"

r = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore', nSP = 20)

r.load_inputfiles(True, False, False)

SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]

rmeteo = r.meteoricR(SPs, 'ms', 1, 4, ret = True)

outpath = "./Data/Calcolo_ricarica_totale"
coordpath = "./Data/Calcolo_ricarica_totale/coord.csv"
geodf = r.georef('recharge', 'rmeteo', coordpath,
         crs = 'epsg:3003', outpath = outpath, setindex = False)

import time
start = time.time()
geodf.to_file(f'{outpath}/prova.shp', index = False)
end = time.time()
print(f'time: {round(end-start, 2)}')

r.georef('recharge', 'rmeteo', coordpath, crs = 'epsg:3003', outpath = outpath)

r.georef('recharge', 'rmeteo', coordpath, crs = 'epsg:3003', outpath = outpath,
         fname = 'rmeteo.geojson', driver = 'GeoJSON')


#con irr


r.load_inputfiles(True, True, False)


coeffs = {
    'E': 0.3,  #Irrigation technique efficiency
    'R': 0.05, #Residual runoff
    'RISP': 1, #1 - fraction of water saved by a change of irrigation technique
    'P': 1     #Percentage of the cell covered by the irrigation
    }
#Path to the input file related to the "special" irrigation district
#Leave it as 'none' if you don't have one
spath = f'{inputpath}/rirrigua_speciale.csv'

r.irrigationR(coeffs, spath)

r.totalR()
r.georef('recharge', 'rmeteo', coordpath, crs = 'epsg:3003', outpath = outpath)
