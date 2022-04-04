# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 16:17:37 2022

@author: paolo
"""

# %% Setup

import os
import pandas as pd
os.chdir('C:/E-OBS-SWB2')

# %% Call the class

from Python.EOBSobject import EOBSobject

# %% Class usage

# 1. Initialize the class

# Define the variables
#Path to the folder where to store the results: outpath
#   it can be a path to a custom folder or direct path to the model folder,
#   as the example below
#   outpath is set as inpath if the field is left untouched
outpath = './Model/swb2_MODELMI/climate_ncfile'
#Path to the folder where the E-OBS data are stored: inpath
inpath = './Data/E-OBS'
#E-OBS code (lowercase) of the variable needed: var
#   All the codes are provided at:
#   https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php
var = 'rr' #daily precipitation sum

f = EOBSobject(inpath, var, outpath)

#If you want to provide directly the path to a single file, set folder to False
# inpath = './Data/E-OBS/file.nc'
# f = EOBSobject(inpath, var, folder = False)

#If the output needs to be used for SWB2
f = EOBSobject(inpath, var, outpath, swb2 = True)

# Load the netcdf file
f.load()

# 2. Cut in space

coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)

# 3. Cut in time

start = 2014
end = 2018

# 4. Cut in space and time

# 5. Generate ArcGRID files
#save_arcgrid

# 6. Perform the operation on multiple E-OBS .nc files

var = ['rr', 'tn', 'tx']
for v in var:
    f = EOBSobject(inpath, v, outpath)
    f.load()
    f.cut_spacetime(coord, start, end)
    f.close()

#The name of the output netcdf file will be composed as var_EOBS_method_year