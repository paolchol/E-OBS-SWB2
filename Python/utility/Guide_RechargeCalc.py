# -*- coding: utf-8 -*-
"""
Guide to the use of class "RechargeCalc"

@author: paolo
"""

# %% Setup

import os
import numpy as np
# import pandas as pd
# import glob

os.chdir('C:/E-OBS-SWB2')

# %% Call the class

from Python.RechargeCalc import RechargeCalc

# %% Class usage

#1. Initialize the class

#Define the variables
startyear = 2014
endyear = 2018
cell_area = 100*100 #m2
#Path to the SWB2 output
swb2path = "./Data/SWB2_output/1Speranza_netinfiltration.nc"
#Path to the input .csv files folder
inputpath = "./Data/Calcolo_ricarica_totale"

r = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore')
r.loadinputfiles()
#if no irrigation or urban recharges are needed, set urb or irr to False
#r.loadinputfiles(urb = False)

#Correct the indicators provided
r.input['ind'].loc[r.input['ind']['distretto'] == 'Muzza', 'distretto'] = 'MUZZA'

#2. Create meteoric recharge dataframe

#Define the stress periods
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

r.meteoricR(SPs)

#3. Create irrigation recharge dataframe

#Define the coefficients needed
coeffs = {
    '1': 0.5,
    '2': 0.7,
    '3': 0.5    #no need of the third if no "special" code 2 is present
    }
#Percentage of irrigation in each stress period
I1 = 0  #%
I2 = 62 #%
I3 = 88 #%
I4 = 0  #%
Is = [I1, I2, I3, I4]
#Path to the input file related to the "special" irrigation district
#Leave it as 'none' if you don't have one
spath = f'{inputpath}/rirrigua_speciale.csv'

r.irrigationR(Is, coeffs, spath)

#4. Create urban recharge dataframe

coeff_urb = 0.15

r.urbanR(coeff_urb)

#5. Create total recharge dataframe

#Launch after computing all the partial recharges
r.totalR()

#Launch directly
#needs the parameters defined before, inside a dictionary or a dataframe
meteopar = {
    'SPs': SPs
    }
irrpar = {
    'coeffs': coeffs,
    'Is': Is,
    'spath': spath
    }
urbpar = {
    'coeff_urb': coeff_urb
    }

r.totalR(meteopar, irrpar, urbpar)

#6. Export

#as .csv
outpath = "./Data/Calcolo_ricarica_totale"
r.export('recharge', 'rtot', outpath = outpath)
#as .shp
#proj: crs in the PROJ4 format. In this case, Monte Mario EPSG:3003
r.georef('recharge', 'rtot', "./Data/Calcolo_ricarica_totale/coord.csv",
         proj = '+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +units=m +no_defs ',
         outpath = outpath)


