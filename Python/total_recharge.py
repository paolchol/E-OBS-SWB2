# -*- coding: utf-8 -*-
"""
Compute the total recharge

- Load the meteorological recharge generated from SWB2
- Divide the meteorological recharge in stress periods (SP). Sum the recharge
    inside the SP and save the sums in a 3D variable: rmeteo(i, nrow, ncol),
    where i is the SP number
- Load databases on land cover, irrigation recharge, municipalities of the
    area and urban recharge for each municipality
- Compute the total recharge summing the databases

@author: paolo
"""

# %% Modules

import os
# import netCDF4 as nc
import numpy as np
import pandas as pd
import glob

# %% Custom functions and classes

os.chdir('C:/E-OBS-SWB2')

#The directory has to be set in ./E-OBS-SWB2 for this to work
from Python.SWB2output import SWB2output
from Python.custom_functions import repeat_list

# %% Setup

#Variables needed
startyear = 2014
endyear = 2018

#Paths
#1. Path to the SWB2 output
swb2path = "./Data/SWB2_output/1Speranza_netinfiltration.nc"
#2. Path to the input .csv files folder
inputpath = "./Data/Calcolo_ricarica_totale"

# %% 0. Load input files

fls = glob.glob(f'{inputpath}/*.csv')
names = ['indicatori', 'ricarica_irrigua', 'ricarica_urbana', 'rirrigua_speciale']
#Assign the paths to the names
k = []
for name in names:
    for i in range(len(fls)):
        if(fls[i].find(name) != -1):
            k += [i]
#Load the files into the variables
ind_df = pd.read_csv(fls[k[0]])
in_rirr = pd.read_csv(fls[k[1]])
in_rurb = pd.read_csv(fls[k[2]])
sp_rirr = pd.read_csv(fls[k[3]])

#Correct the indexes in ind_df
ind_df.loc[ind_df['distretto'] == 'Muzza', 'distretto'] = 'MUZZA'

# %% 1. Meteoric recharge

#Stress period definition
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

#Set up the class
f = SWB2output(swb2path)
#Check the units (should be in inches)
f.metadata['units'] #inches
#Return the SP sum directly in m/s
rmeteo = f.SP_sum(SPs, units = 'ms')
#Close the netCDF file to save memory
f.close()

# %% 2. Irrigation recharge

#Percentage of irrigation in each stress period
I1 = 0  #%
I2 = 62 #%
I3 = 88 #%
I4 = 0  #%
#Create the list and repeat it n times, where n is the number of years
Is = [I1, I2, I3, I4]
Is = repeat_list(Is, endyear-startyear+1, True)

#Define the coefficients needed to calculate the irrigation recharge
coeff1 = 0.5
coeff2 = 0.7
coeff3 = 0.5 #Defined for 'distretto' with code 'speciale' = 2
cell_area = 100*100 #m2

#Calculate the irrigated area
in_rirr.insert(len(in_rirr.columns),'area', 0)
for i, distr in enumerate(in_rirr['distretto'], 0):
    area = sum((ind_df['distretto'] == distr) & (ind_df['zona_agricola'] == 1))*cell_area
    in_rirr.loc[i, 'area'] = area
del area
#Calculate the discharge Q in m/s
in_rirr.insert(len(in_rirr.columns),'Q_ms', 0)
in_rirr['Q_ms'] = in_rirr['portata_concessa_ls'] * 0.001 / in_rirr['area']

#Calculate the dataframe of irrigation recharge
#rirr will be in m/s
rirr = ind_df.loc[:,('indicatore', 'distretto', 'zona_agricola')]
for i, I in enumerate(Is, 0):
  for j, distr in enumerate(in_rirr['distretto'], 0):
    sp = in_rirr.loc[j, 'speciale'] #'special' code
    if(sp != 1):
        Q = in_rirr['Q_ms'][j]
        cond = (rirr.loc['distretto'] == distr) & (rirr['zona_agricola'] == 1)
        if f'SP{i+1}' not in rirr.columns:
            rirr.insert(len(rirr.columns), f'SP{i+1}', 0)        
        rirr.loc[cond, f'SP{i+1}'] = Q * I * coeff1 * coeff2 if sp != 2 else Q * I * coeff1 * coeff2 * coeff3
#Assign the provided "special" recharge to the "special" districts
special = in_rirr.loc[in_rirr['speciale'] == 1, 'distretto']
for s in special:
    cond = (rirr['distretto'] == s) & (rirr['zona_agricola'] == 1)
    rirrcol = rirr.columns[3:]
    spcol = sp_rirr.columns[1:]
    rirr.loc[cond, rirrcol] = sp_rirr.loc[sp_rirr['distretto'] == s, spcol].values

#To clean up unnecessary values after this section
del Q, cond, sp, i, j, distr, spcol, rirrcol

# %% 3. Urban recharge


# %% 4. Total recharge


# %% 5. Export the results

#Total recharge

#Partial recharges
# - Meteoric recharge

# - Irrigation recharge

# - Urban recharge
