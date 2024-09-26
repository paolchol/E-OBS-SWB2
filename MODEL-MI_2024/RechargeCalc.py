'''
Calculate total recharge
Based on Guide_RechargeCalc.py by @paolochol

Input files:
- net_infiltration.nc (output from nc_outputSWB2.py)
- indicatori_update.csv
- NEWextractions.csv
- NEWricarica_irrigua.csv
- NEWirrigua_speciale.csv

Class(es):
- EOBS-SWB2/Python/SWB2output.py
- EOBS-SWB2/Python/RechargeCalc.py

'''

import os
import numpy as np

main_dir = 'C:/Users/user/OneDrive - Politecnico di Milano/' 
os.chdir(os.path.join(main_dir,'GitHub/E-OBS-SWB2'))
from Python.RechargeCalc import RechargeCalc

#Define the variables
startyear = 2019
endyear = 2022
cell_area = 100*100 #m2

#Path to the SWB2 output
swb2path = os.path.join(main_dir,'SWB2/MODEL-MI/output/R1_net_infiltration__2019-01-01_2022-12-31__338_by_660.nc')
#Path to the input .csv files folder
inputpath = os.path.join(main_dir,'GitHub/E-OBS-SWB2/MODEL-MI_2024/Data_modified')

#1. Initialize class       
r = RechargeCalc(startyear, endyear, cell_area, uniqueid = 'indicatore', nSP = 16)

#Load the input files needed
r.load_inputfiles(swb2path, inputpath)
#Access input files via: r.input['ind']

#2. Create meteoric recharge dataframe

#Define the stress periods lengths
SP1 = 90   #days, 01/01 - 31/03
SP2 = 76   #days, 01/04 - 15/06
SP3 = 92   #days, 16/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]

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
spath = f'{inputpath}/rirrigua_speciale_new.csv'

r.irrigationR(coeffs, spath)

'''
#Multiple coefficients --> APPLY?
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
'''

#4. Create urban recharge dataframe
#Compute the urban recharge as a fraction of the pumped volumes. It is due to the losses from the extraction pumps and pipes.
coeff_urb = 0.15  #Coefficient to apply to the extraction

#Specify the condition you want to apply by providing the columns you want to consider (col) and which values they need (valcol).
#Specify also which operation you want to check between the conditions

# (indicatori['land_cover'] == 123) OR (indicatori['land_cover'] == 124) AND (indicatori['zona_urbana'] == 1)
col = ['land_cover', 'land_cover', 'zona_urbana']
valcol = [123, 124, 1]
option = [0, 1] #0: OR, 1: AND

#Call urbanR()
r.urbanR(coeff_urb, col = col, valcol = valcol, option = option)

'''
#If you want to return the areas of the cells in the specified condition
#you can add areas = True to the function
r.urbanR(coeff_urb, col = col, valcol = valcol, option = option, areas = True)
urb = r.get_df('input', 'urb')
'''

#5. Create total recharge dataframe
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

#Export as .csv with coords
outpath = os.path.join(main_dir,'GitHub/E-OBS-SWB2/MODEL-MI_2024/Output')
r.export('recharge', 'rtot', outpath = outpath, outname = 'Prova1',
         withcoord = True, coordpath = f'{inputpath}/coord.csv')

'''
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
'''