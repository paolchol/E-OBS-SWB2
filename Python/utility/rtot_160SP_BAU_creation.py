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
inputpath = "./Stefano"

r = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore', nSP = 160)
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
spath = f'{inputpath}/rirrigua_speciale.csv'

r.irrigationR(coeffs, spath)

r.totalR()

r.georef('recharge', 'rtot', "./Stefano/coord.csv",
         proj = '+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +units=m +no_defs ',
         outpath = "./Stefano", outname = 'rtot_160SP_BAU', dropcoord = True)
#35 minutes needed
