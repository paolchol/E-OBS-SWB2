# -*- coding: utf-8 -*-
"""
Climate change scenarios creation for SWB2

- Load the raster with the number of the station as a dataframe
    (be careful on using the same numbers)
- Load the climate change scenarios
- For each climate change scenario:
    For each station:
        - Open the station dataframe
        - Get the ith day
        - Assign the value to the corresponding tile of the TIFF

@author: paolo
"""

import glob
import numpy as np
import os
import pandas as pd
from PIL import Image

#Set the working directory
os.chdir("C:/E-OBS-SWB2")

from Python.custom_functions import save_ArcGRID

#Set all the variables needed
xll = 1
yll = 1
cellsize = 100




#Load the GeoTIFF and set it as a numpy.array
im = Image.open('./Data/Raster/ModelMI_ID_stazioni.tif')

#Set up the folders in which the climate change scenarios are stored
#Example:
# CCsc1 > station1.csv, station2.csv
# CCsc2 > station1.csv, station2.csv
cclist = [x[0] for x in os.walk('./Data/CCscenarios/')]

start = 2014
end = 2018

for ccpath in cclist:
    stations = glob.glob(f'{ccpath}/*.csv')
    for day in range(start, end): #get the row
        imarray = np.array(im)
        for stpath in stations:
            stcode = 1 #get this from the filename and by interrogating a database of station name+code
            station = pd.read_csv(stpath)
            imarray[imarray == stcode] = station.loc[day, 'value'] #insert the value in the raster
        fname = ''
        save_ArcGRID(imarray, fname)
        
    

