# -*- coding: utf-8 -*-
"""
SW2 output handling
Creation of a class to handle the SWB2 output

@author: paolo
"""

# %% Setup

import os
import netCDF4 as nc
import numpy as np
import pandas as pd

os.chdir('C:/E-OBS-SWB2')
from Python.custom_functions import save_ArcGRID

# swb2path = "./Model/swb2_MODELMI/output/ModelMI_net_infiltration__2014-01-01_2018-12-31__338_by_660.nc"
swb2path = "./Data/SWB2_output/p1_ModelMI_net_infiltration__2014-01-01_2018-12-31__338_by_660.nc"

# %% Class definition

class SWB2output():
    
    #Custom functions needed
    def getkeys(dict):
        #got from: https://www.geeksforgeeks.org/python-get-dictionary-keys-as-a-list/
        # return list(dict.keys()) other method
        return [*dict]
    
    def getdates(string, n):
        from datetime import date, timedelta
        datestr = string.split(' ')[2]
        y, m, d = datestr.split('-')
        start = date(int(y), int(m), int(d))
        end = start + timedelta(days = n-1)
        return start, end
    
    #Actual class functions definition
    def __init__(self, path):
        #Load the output netCDF file
        #Extract its main metadata
        start, end = getdates(f.netCDF.variables['time'].units,
                              f.netCDF.variables['time'].shape[0])
        self.path = path
        self.netCDF = nc.Dataset(path)
        self.metadata = {
            "start_date": start,
            "end_date": end,
            "variables": getkeys(self.netCDF.variables),
            "main_variable": getkeys(self.netCDF.variables)[3], 
            "units": f.netCDF[getkeys(f.netCDF.variables)[3]].units
            }
    
    def sumtot(self, outpath = 'None', name = 'name'):
        #Returns the sum over the whole time period
        #If a path is provided as outpath, an ArcASCII GRID file is produced
        variable = self.metadata["main_variable"]
        self.sumtotdf = np.sum(np.ma.getdata(self.netCDF[variable][:,:,:]), axis = 0)*0.0254 #meters
        
        if(outpath != 'None'):
            df = pd.DataFrame(self.sumtotdf)
            size = round(np.ma.getdata(self.netCDF['x'][1]).item() - np.ma.getdata(self.netCDF['x'][0]).item())
            xll = round(np.ma.getdata(self.netCDF['x'][0]).item()) - size/2
            yll = round(np.ma.getdata(self.netCDF['y'][-1]).item()) - size/2
            #name: prendere ModelMI. nel caso in cui sia p1_ prenderà p1
            fname = f"{outpath}/{name}_{variable}_sum_tot.asc"
            save_ArcGRID(df, fname, xll, yll, size, -9999)
            print(f'ArcGRID saved in {outpath} as: {name}_{variable}_sum_tot.asc')
        return self.sumtotdf
    
    def SP_sum(self, SPs):
        print(SPs)
    
    def close(self):
        self.netCDF.close()

# %% Use the class

#Create the object
f = SWB2output(swb2path)

#Print the metadata
f.metadata #print all the metadata available
print(f.metadata['end_date']) #print one variable of the metadata

#Sum the components
sumtot = f.sumtot('./Data') #to assign the sum to a new variable
f.sumtotdf #in any case, the sum is saved inside the object, as sumtotdf

#Close the netCDF file
f.close()



#Scrivere come utilizzare la classe man mano

# 		def metadata(self):
# 			#returns the output metadata, name of the variable ecc
# 		
#         def sum(self, varname):
# 			#returns the sum of the variable specified over the years
# 			#could varname be exracted from the netCDF metadata in some way?
# 			#print the name of the variable summed
# 		def plot(self):
# 			#plot the variable' sum. As a raster
# 		def histogram(self):
# 			#plot a histogram of the variable sum
# 		def SP_sum(self, SP_bounds):
# 		#divide the dataset into stress periods (SP) and create the sum

# %% How to get information from a netCDF file

#Variables is a dictionary
#Get the names of the variables
var = f.netCDF.variables #prints all the variables and their metadata
keys = f.netCDF.variables.keys()
print(keys)

def getkeys(dict):
    #got from: https://www.geeksforgeeks.org/python-get-dictionary-keys-as-a-list/
    # return list(dict.keys()) other method
    return [*dict]

keys = getkeys(f.netCDF.variables)
getkeys(f.netCDF.variables)[3] #only the main variable

#Start and end dates

f.netCDF.variables['time'].units #starting date
f.netCDF.variables['time'].shape[0] #time length

def getdates(string, n):
    from datetime import date, timedelta
    datestr = string.split(' ')[2]
    y, m, d = datestr.split('-')
    start = date(int(y), int(m), int(d))
    end = start + timedelta(days = n-1)
    return start, end

string = f.netCDF.variables['time'].units
n = f.netCDF.variables['time'].shape[0]
start, end = getdates(string, n)

#Units
f.netCDF[getkeys(f.netCDF.variables)[3]].units



