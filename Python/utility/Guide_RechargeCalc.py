
# -*- coding: utf-8 -*-
"""
Guide to the use of class "RechargeCalc"

@author: paolo
"""

# %% Setup

import os
import numpy as np

os.chdir('C:/E-OBS-SWB2')

# %% Call the class

from Python.RechargeCalc import RechargeCalc

# %% Class usage

#0. Set up your input files

#1. Initialize the class

#Define the variables
startyear = 2014
endyear = 2018
cell_area = 100*100 #m2
#Path to the SWB2 output
swb2path = "./Data/SWB2_output/VersioneFINALE_net_infiltration.nc"
#Path to the input .csv files folder
inputpath = "./Data/Calcolo_ricarica_totale"

r = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore', nSP = 20)

#Load the input files needed
r.load_inputfiles()
#if no irrigation or urban recharges are needed, set meteo, urb or irr to False
r.load_inputfiles(urb = False)

#Once you loaded the input files, you can access them via
r.input['ind']
r.input['irr']
#You can then perform any operation you would on dataframes,
#for example correct a wrong value provided
r.input['ind'].loc[r.input['ind']['distretto'] == 'Muzza', 'distretto'] = 'MUZZA'

#2. Create meteoric recharge dataframe

#Define the stress periods
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]

#Launch meteoricR() by providing the SPs, the units wanted and the starting row and column
#of the "indicatori" file to correctly associate the cells from the SWB2 output to the 
#cells of the desired area
r.meteoricR(SPs, 'ms', 1, 4)

#3. Create irrigation recharge dataframe

#Define the coefficients needed
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

#Multiple coefficients
#Apply a different set of coefficients to a selected list of stress periods. In
#this example, E = 0, R = 0, RISP = 0 and P = 1 are applied to the stress periods
#SP3 and SP4
multicoeffs = {
    'E': [0.3, 0],  #Irrigation technique efficiency
    'R': [0.05, 0], #Residual runoff
    'RISP': [1, 0], #1 - fraction of water saved by a change of irrigation technique
    'P': [1, 1]     #Percentage of the cell covered by the irrigation
    }
splist = ['SP3', 'SP4']

r.irrigationR(multicoeffs, spath, multicoeff = True, splist = splist)

#4. Create urban recharge dataframe

#The class gives the possibility to compute the urban recharge as a fraction
#of the pumped volumes. It is due to the losses from the extraction pumps and pipes.
#To generate the urban recharge dataframe you can set up custom conditions
#to identify in which cells you want to compute the urban recharge.
#You have to specify the coefficient to apply to the extraction
coeff_urb = 0.15
#Then, you have to specify the condition you want to apply by providing the
#columns you want to consider (col) and which values they need (valcol). In options
#you have to specify which operation you want to check between the conditions
#0 if OR (| in python), 1 if AND (& in python). col and valcol need to be lists
#even if they contain a single value, while option it's not needed in this case
#The example below means:
# (indicatori['land_cover'] == 123) OR (indicatori['land_cover'] == 123) AND (indicatori['zona_urbana'] == 1)
col = ['land_cover', 'land_cover', 'zona_urbana']
valcol = [123, 124, 1]
option = [0, 1] #0: OR, 1: AND

#Call urbanR()
r.urbanR(coeff_urb, col = col, valcol = valcol, option = option)

#If you want to return the areas of the cells in the specified condition
#you can add areas = True to the function
r.urbanR(coeff_urb, col = col, valcol = valcol, option = option, areas = True)
urb = r.get_df('input', 'urb')

#5. Create total recharge dataframe

#Launch after computing all the needed partial recharges
r.totalR()

#Launch directly without computing the partial recharges one by one.
#needs the parameters defined before, inside a dictionary or a dataframe
meteopar = {
    'SPs': SPs
    }
irrpar = {
    'coeffs': coeffs,
    'spath': spath
    }
urbpar = {
    'coeff': coeff_urb,
    'col': col,
    'valcol': valcol,
    'option': option
    }

r.totalR(meteopar, irrpar, urbpar)

#6. Export

#Export as .csv
outpath = "./Data/Calcolo_ricarica_totale"
r.export('recharge', 'rtot', outpath = outpath)
#You can include the X and Y coordinates in the CSV file
r.export('recharge', 'rtot', outpath = outpath, outname = 'prova',
         withcoord = True, coordpath = "./Data/Calcolo_ricarica_totale/coord.csv")

#Export as .shp
#proj: crs in the PROJ4 format. In this case, Monte Mario EPSG:3003
r.georef('recharge', 'rtot', "./Data/Calcolo_ricarica_totale/coord.csv",
         crs = 'epsg:3003',
         # proj = '+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +units=m +no_defs ',
         outpath = outpath)
#If you don't want to keep the X and Y columns in the exported dataframe,
#set dropcoord as True
r.georef('recharge', 'rtot', "./Data/Calcolo_ricarica_totale/coord.csv",
         crs = 'epsg:3003',
         outpath = outpath, outname = 'rtot_noXY', dropcoord = True)

#7. Extract or modify values

#Extract
#Use get_df and .copy() to create an exact copy of the DataFrame. In this 
# way any modification you may do on extracted will not affect the DataFrame
# stored in the RechargeCalc object
extracted = r.get_df('recharge', 'rtot').copy()

#Modify - Examples
#Modify the urban cells values multiplyng them by 2
r.modify_recharge('recharge', 'rtot', 2,
                  col = 'zona_urbana', valcol = 1)

#Modify the urban cells inside a specified municipality, multiplying them by 3
r.modify_recharge('recharge', 'rtot', 3,
                  single_cond = False, multi_cond = True,
                  col = ['zona_urbana', 'nome_com'], valcol = [1, 'MILANO'])

#Modify the urban cells inside a specified municipality and for a specified
# land cover, multiplying them by 4
r.modify_recharge('recharge', 'rtot', 4,
                  single_cond = False, multi_cond = True,
                  col = ['zona_urbana', 'nome_com', 'land_cover'], valcol = [1, 'MILANO', 121])

#You can then export the desired modified dataframe in the same ways as
# explained in Section 6
outpath = "./Stefano"
r.georef('recharge', 'rtot', "./Data/Calcolo_ricarica_totale/coord.csv",
         crs = 'epsg:3003',
         outpath = outpath, outname = 'rtot_2urb', dropcoord = True)
#All the modifies add up. To return to the original recharge, re-run the 
#function used to obtain the value in the first place. In this case r.totalR()
r.totalR()
