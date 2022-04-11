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
outpath = './Data'
#Path to the folder where the E-OBS data are stored: inpath
inpath = './Data/E-OBS'
#E-OBS code (lowercase) of the variable needed: var
#   All the variable codes are provided at:
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

# 2. Generate netcdf files
#The name of the output netcdf file will be composed as var_EOBS_method_year

# 2.1 Cut in space

coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)
f.cut_space(coord)

#To add more cells
f.set_fname(f'{outpath}/rr_morecells.nc')
#just to provide a custom name to distinguish the files, not needed for the code to work
f.cut_space(coord, contourcell = 2)

# 2.2 Cut in time

start = 2014
end = 2018
f.cut_time(start, end)

#The default option will generate one file for each year. If you want to generate
#a single file, set option as 'bundle'
f.cut_time(start, end, option = 'bundle')

#With the bundle option, you can also cut between given days
#You need to set day as True and provide start and end as datetime.date objects
from datetime import date
start = date(2011, 12, 30)
end = date(2017, 7, 15)
f.cut_time(start, end, option = 'bundle', day = True)

# 2.3 Cut in space and time

#All options available for cut_space and cut_time are also available for cut_spacetime
f.cut_spacetime(coord, start, end)

# 2.4 Keep the raw file

#You can also save the file as it is
#This will only change the metadata or other things as the name of the main variable

f.save_netcdf(method = 'raw')

# 3. Generate daily ArcGRID files

#Specify the method: 'cut_space', 'cut_time', 'cut_spacetime'
#Provide the necessary information: coord, start, end

f.save_arcgrid()

# 4. Perform the operation on multiple E-OBS .nc files

var = ['rr', 'tn', 'tx']
#Also outpath and inpath could be provided in lists
for v in var:
    f = EOBSobject(inpath, v, outpath)
    f.load()
    f.cut_spacetime(coord, start, end)
    f.close()

