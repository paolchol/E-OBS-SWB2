# EOBSobject class guide

## Introduction

E-OBS data are bla bla...
Link to the data
Explain netcdf format
Link to the netcdf project page

## Code

### Setup

Import the necessary modules and set up the working directory
'''python
import os
import pandas as pd
os.chdir('C:/E-OBS-SWB2')
'''

Import the class
'''python
from Python.EOBSobject import EOBSobject
'''
'Python.' is only needed if the EOBSobject.py file is in another folder as in this repository. If the file is in the same folder as the main only 'from EOBSobject' is needed.

### 1. Initialize the class

#### Define the variables
Path to the folder where the E-OBS data are stored: 'inpath'. *This is a required parameter*
Path to the folder where to store the results: 'outpath'
E-OBS code (lowercase) of the variable needed: 'var' *This is a required parameter*
'outpath' can be a path to a custom folder or direct path to the model folder, as the example below.
'outpath' is set as 'inpath' if the field is left untouched.
'''python
inpath = './Data/E-OBS'
outpath = './Model/swb2_MODELMI/climate_ncfile'
var = 'rr' #daily precipitation sum
'''
All the E-OBS variable codes are provided at: https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php
The files in the inpath folder which are not the one you want to load shouldn't have the name of the variable (var) at the start of their file name. The code could identify them as the E-OBS file and thus not work. Example: *rr_otherfile.pdf* shouldn't be in the folder.

#### Create the object
'''python
f = EOBSobject(inpath, var, outpath)
'''

If you want to provide directly the path to a single file, set 'folder' to 'False' when creating the object.
'''python
inpath = './Data/E-OBS/file.nc'
f = EOBSobject(inpath, var, folder = False)
'''

The 'swb2' option has to be set to True if the output needs to be used for SWB2.
'''python
f = EOBSobject(inpath, var, outpath, swb2 = True)
'''

#### Load the netcdf file
'''python
f.load()
'''