# -*- coding: utf-8 -*-
"""
Guide of the "SWB2output" class for SWB2 output handling
 In section 'Guide to the class usage', instructions are provided 

@author: paolo
"""

# %% Setup

import os
import netCDF4 as nc
import numpy as np
import pandas as pd

os.chdir('C:/E-OBS-SWB2')

# %% Call the class

from Python.SWB2output import SWB2output

# %% Guide for the class use

#Path to an example netCDF file produced by SWB2
swb2path = "./Data/SWB2_output/p1_ModelMI_net_infiltration__2014-01-01_2018-12-31__338_by_660.nc"

#1. Create the object
f = SWB2output(swb2path)

#2. Print the metadata
f.metadata #print all the metadata available
print(f.metadata['end_date']) #print one variable of the metadata

#3. Sum the main variable over the whole period
#The sum can be assigned to a new variable
sumtot = f.sumtot('./Data')
f.sumtotdf
#In any case, the sum is saved inside the object, as .sumtotdf

#4. Sum over the stress periods
#Stress period definition
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 
SP4 = 107  #days, - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

#The sum can be assigned to a new variable
SPsum = f.SP_sum(SPs)
#Specify outpath if you want to save the ArcASCII GRID files
f.SP_sum(SPs, outpath = 'some path')

#5. Close the netCDF file
#Important to save memory
f.close()

