#questo script verrà poi messo al posto del main che c'è ora, nella cartella principale E-OBS-SWB2

#What this main is and what it does

#Explanation on how to use the main

## 0. Setup

import os

#Set the directory as the position of the script
#os.dir something

#Launch custom functions and needed modules (as source in R)

#Create the basic folders: Data, Export

## 1. SWB2 input creation

### 1.1 General setup

#Write the path to the existing folder where E-OBS data are
#Default path: "./Data/E-OBS"
eobs_path = "./Data/E-OBS"

#Which type of input do you need?
# - write 'ASCII' for ArcASCII GRID
# - write 'netCDF' for netCDF
select_input = 'netCDF' 

#In the script inputgen, place an if statement checking select_input

#Where do you want the output to be pruced?
#Write the path to the existing folder where you want it to be produced
#Default path: "./Export/netCDF"
export_path = "./Export/netCDF"

#In the script inputgen, search for the folder, if it doesn't exist, create it
#If it can't be created, return a message saying that the folder is missing

### 1.2 Launch the input generation

#Specify the desired projection of the export
#If no reprojection is needed, set reproj as False
#Valid only for netCDF files, ArcASCII export will be kept in WGS84 reference system
reproj = True
utmzoneN = 32
utmzoneL = 'N'

#in inputgen script: set an if statement for reproj

#Set the extension of the final area
#Provide the path to a shapefile, a csv, a txt or a pandas DataFrame containing the coordinates of the
# maximum extension points of the area over which "cut" the E-OBS file
extension = './Data/ancillary_grid_extentions_SWB.txt'

#in the inputgen script, check the extention of the file, and based on that select the operation
#needed to load the file and extract maxlat, maxlon ecc.

#Launch the script
#equivalent of source('./packages/input_gen.py')

## 2. SWB2 output handling

#Launch the script
#equivalent of source('./packages/input_gen.py')