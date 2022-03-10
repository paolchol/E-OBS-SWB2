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

#Paths
#1. Path to the SWB2 output
swb2path = "./Data/SWB2_output/1Speranza_netinfiltration.nc"

# %% Obtain the stress periods sums

#Stress period definition
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 
SP4 = 107  #days, - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

f = SWB2output(swb2path)
SPsum = f.SP_sum(SPs)
f.close()

# %% Load and organize the needed datasets



# %% Compute the total recharge


