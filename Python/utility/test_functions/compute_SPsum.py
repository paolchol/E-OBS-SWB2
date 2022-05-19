# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 11:50:47 2022

@author: paolo
"""

import os
# import netCDF4 as nc
import numpy as np
import pandas as pd
import glob

# %% Custom functions and classes

os.chdir('C:/E-OBS-SWB2')

#The directory has to be set in ./E-OBS-SWB2 for this to work
from Python.SWB2output import SWB2output

# %% Sum over the stress periods

path = "./Data/SWB2_output/0Impervious_net_infiltration.nc"
f = SWB2output(path)

SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 
SP4 = 107  #days, - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)

_ = f.SP_sum(SPs, outpath = r".\Export\ASCII\RMeteo_SP_inch\impervious0")