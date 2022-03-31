# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 16:17:37 2022

@author: paolo
"""

# %% Setup

import os
os.chdir('C:/E-OBS-SWB2')

# %% Call the class

from Python.EOBSobject import EOBSobject

# %% Class usage

# 1. Initialize the class

# Define the variables
#Path to the folder where to store the results: outpath
#   it can be a path to a custom folder or direct path to the model folder,
#   as the example below
outpath = './Model/swb2_MODELMI/climate_ncfile'
#Path to the folder where the E-OBS data are stored: inpath
inpath = './Data/E-OBS'
#E-OBS code (lowercase) of the variable needed: var
#   All the codes are provided at:
#   https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php
var = 'rr' #daily precipitation sum

f = EOBSobject(inpath, outpath, var)
f.load() #load the netcdf file

#If you want to provide directly the path to a single file, set folder to False
# inpath = './Data/E-OBS/file.nc'
# f = EOBSobject(inpath, var, folder = False)

#outpath is set as inpath if the field is left untouched


