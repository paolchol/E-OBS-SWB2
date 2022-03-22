# -*- coding: utf-8 -*-
"""
Guide to the use of class "RechargeCalc"

@author: paolo
"""

import os
# import netCDF4 as nc
import numpy as np
import pandas as pd
import glob

os.chdir('C:/E-OBS-SWB2')
from Python.RechargeCalc import RechargeCalc

#Define the variables

startyear = 2014
endyear = 2018
cell_area = 100*100 #m2
#Path to the SWB2 output
swb2path = "./Data/SWB2_output/1Speranza_netinfiltration.nc"
#Path to the input .csv files folder
inputpath = "./Data/Calcolo_ricarica_totale"

# Initialize the class

r = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore')
r.loadinputfiles()

#Correct the indicators provided
r.input['ind'].loc[r.input['ind']['distretto'] == 'Muzza', 'distretto'] = 'MUZZA'

# Create meteoric recharge dataframe

#Define the stress periods
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

r.meteoricR(SPs)

# Create irrigation recharge dataframe

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

# Create urban recharge dataframe

coeff_urb = 0.15

r.urbanR(coeff_urb)

# Create total recharge dataframe

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

# Export


# %% Drafts

getkeys(r.recharges)

r.info['id']
r.recharges['rmeteo'].columns


def findSPcol(col, ind):
    print(ind)
    names = [ind]
    for name in col:
        if name.find('SP') != -1:
            names += [name]
    return names


findSPcol(r.recharges['rmeteo'].columns, r.info['id'])

r.input['ind'].dtypes

cc = r.input['ind'].astype({'indicatore': str})
cc.dtypes

max(r.input['ind']['column'])


