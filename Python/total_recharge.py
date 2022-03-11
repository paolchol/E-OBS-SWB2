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
import netCDF4 as nc
import numpy as np
import pandas as pd

# %% Custom functions and classes

os.chdir('C:/E-OBS-SWB2')

#The directory has to be set in ./E-OBS-SWB2 for this to work
from Python.SWB2output import SWB2output

# %% Setup

#Variables needed
startyear = 2014
endyear = 2018

#Paths
#1. Path to the SWB2 output
swb2path = "./Data/SWB2_output/1Speranza_netinfiltration.nc"

# %% 1. Meteoric recharge

#Stress period definition
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

f = SWB2output(swb2path)
f.metadata['units'] #inches
#Return the SP sum directly in m/s
rmeteo = f.SP_sum(SPs, units = 'ms')
f.close()

# %% 2. Irrigation recharge

#Percentage of irrigation in each stress period
I1 = 0
I2 = 62 #%
I3 = 88 #%
I4 = 0
Is = [I1, I2, I3, I4]

#Calculate the correct PC

#place 'distretto' = 0 as 'null'

#Nuovo metodo senza ciclo, fatto alla fine!

# for SP in Is:
#     for i in rows:
#         for j in columns:
#             ID = int(f'{i}0{j}')
#             cond = ind['distretto'][where 'indicatore' == ID]
#             if cond != 'null':
#                 rirr[SP, i, j] = 
#                 PC['Distretto' == cond, 'portata_concessa']*Is[SP]*0.7*0.5
#                  if ind['Zona_irrigua'] == 1
#                  else 0

 #mettere a punto il ciclo
 #fare un prova e vedere bene quanto ci mette a ottenere rirr completo

# %% 3. Urban recharge


# %% 4. Total recharge


# %% 5. Export the results

#Total recharge

#Partial recharges
# - Meteoric recharge

# - Irrigation recharge

# - Urban recharge
