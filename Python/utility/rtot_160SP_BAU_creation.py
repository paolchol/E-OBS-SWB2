# -*- coding: utf-8 -*-
"""
Created on Wed May  4 10:44:06 2022

@author: paolo
"""

import os
import numpy as np

os.chdir('C:/E-OBS-SWB2')

from Python.RechargeCalc import RechargeCalc

startyear = 2014
endyear = 2018
cell_area = 100*100 #m2
#Path to the SWB2 output
swb2path = "./Data/SWB2_output/VersioneFINALE_net_infiltration.nc"
#Path to the input .csv files folder
inputpath = "./Stefano/360SP"

r = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore', nSP = 360)
r.load_inputfiles(urb = False)

SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]

r.meteoricR(SPs)

coeffs = {
    'E': 0.3,  #Irrigation technique efficiency
    'R': 0.05, #Residual runoff
    'RISP': 1, #1 - fraction of water saved by a change of irrigation technique
    'P': 1     #Percentage of the cell covered by the irrigation
    }
#Path to the input file related to the "special" irrigation district
#Leave it as 'none' if you don't have one
# spath = f'{inputpath}/rirrigua_speciale.csv'

# r.irrigationR(coeffs, spath)
r.irrigationR(coeffs)

r.totalR()

geodf = r.georef('recharge', 'rtot', f'{inputpath}/coord.csv',
                 crs = 'epsg:3003',
                 # proj = '+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +units=m +no_defs'
                 outpath = inputpath, outname = 'rtot_360SP_BAU', dropcoord = True)
#35 minutes needed for 160-sp

geodf.drop(['X', 'Y'], 1, inplace = True)

import geopandas
import fiona
import time

def to_file(df, filename, driver="ESRI Shapefile", schema=None, **kwargs):
    if schema is None:
        schema = geopandas.io.file.infer_schema(df)
    filename = os.path.abspath(os.path.expanduser(filename))
    with fiona.Env():
        with fiona.open(filename, 'w', driver=driver, crs=df.crs,
                        schema=schema, **kwargs) as colxn:
            # small adaptation to original code to split
            # materializing of the iterfeatures generator and fiona's writerecords
            records = list(df.iterfeatures())
            colxn.writerecords(records)


start = time.time()
to_file(geodf, f'{inputpath}/test1.shp')
end = time.time()
print(f'*** With custom to_file: {end-start} s')

geodf.crs

start = time.time()
geodf.to_file(f'{inputpath}/test2.zip', driver = 'ESRI Shapefile')
end = time.time()

print(f'*** With standard geodf.to_file: {end-start} s')

geodf.info()

#si possono salvare più shp e fare un merge successivo su qGIS


#Try to export a csv with x and y coordinates and georeference it in QGIS
#3 hours for 360SP

r.export('recharge', 'rtot', outpath = inputpath,
         withcoord = True, coordpath = f'{inputpath}/coord.csv')

