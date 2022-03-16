# -*- coding: utf-8 -*-
"""
Class "SWB2output" for SWB2 netCDF output handling
It can be called using:
from Python.SWB2output import SWB2output

The working directory has to be set in ./E-OBS-SWB2 for this to work

@author: paolo
"""

import os
import netCDF4 as nc
import numpy as np
import pandas as pd

# os.chdir('C:/E-OBS-SWB2')
from Python.custom_functions import save_ArcGRID
from Python.custom_functions import leap
from Python.custom_functions import getdates
from Python.custom_functions import getkeys

class SWB2output():
    
    def __init__(self, path):
        #Load the output netCDF file
        #Extract its main metadata
        
        self.path = path
        self.netCDF = nc.Dataset(path)
        start, end = getdates(self.netCDF.variables['time'].units,
                              self.netCDF.variables['time'].shape[0])
        self.metadata = {
            "start_date": start,
            "end_date": end,
            "variables": getkeys(self.netCDF.variables),
            "main_variable": getkeys(self.netCDF.variables)[3], 
            "units": self.netCDF[getkeys(self.netCDF.variables)[3]].units,
            "nodata_value": self.netCDF[getkeys(self.netCDF.variables)[3]]._FillValue
            }
    
    def print_md(self):
        #Prints the metadata
        print(self.netCDF)
    
    def extract(self):
        #Returns the main variable of the output as a numpy array
        variable = self.metadata['main_variable']
        return np.ma.getdata(self.netCDF[variable][:, :, :])

    def sumtot(self, outpath = 'none', name = 'name'):
        #Returns the sum over the whole time period
        #If a path is provided as outpath, an ArcASCII GRID file is produced
        variable = self.metadata["main_variable"]
        self.sumtotdf = np.sum(np.ma.getdata(self.netCDF[variable][:,:,:]), axis = 0)*0.0254 #meters
        
        if(outpath != 'none'):
            df = pd.DataFrame(self.sumtotdf)
            size = round(np.ma.getdata(self.netCDF['x'][1]).item() - np.ma.getdata(self.netCDF['x'][0]).item())
            xll = round(np.ma.getdata(self.netCDF['x'][0]).item()) - size/2
            yll = round(np.ma.getdata(self.netCDF['y'][-1]).item()) - size/2
            #name: prendere ModelMI. nel caso in cui sia p1_ prenderà p1
            fname = f"{outpath}/{name}_{variable}_sum_tot.asc"
            save_ArcGRID(df, fname, xll, yll, size, -9999)
            print(f'ArcGRID saved in {outpath} as: {name}_{variable}_sum_tot.asc')
        return self.sumtotdf
    
    def SP_sum(self, SPs, outpath = 'none', name = 'name', units = 'none'):
        #Perform a sum of the main variable over the stress periods provided
        # as SPs
        #Returns a 3D variable containing all the sums as the 0 index
        #If outpath is specified, saves an ArcASCII GRID file for each sum
        
        variable = self.metadata['main_variable']
        units = self.metadata['units'] if units == 'none' else units
        print(f'Performing the sum of {variable} over the stress periods provided')
        print(f'Output unit measure: {units}')
        
        #Get the necessary variables
        starty = self.metadata['start_date'].year
        endy = self.metadata['end_date'].year
        if(outpath != 'none'):
            size = round(np.ma.getdata(self.netCDF['x'][1]).item() - np.ma.getdata(self.netCDF['x'][0]).item())
            xll = round(np.ma.getdata(self.netCDF['x'][0]).item()) - size/2
            yll = round(np.ma.getdata(self.netCDF['y'][-1]).item()) - size/2
        
        period = range(starty, endy+1)
        s, e, k = 0, 0, 0
        #Create the 3D variable
        var3d = np.zeros((len(period)*len(SPs), self.netCDF[variable].shape[1], self.netCDF[variable].shape[2]))
        
        for y in period: #Extract a single year
            e += leap(y)
            year = np.ma.getdata(self.netCDF[variable][s:e, :, :])
            base = 0
            for i, SP in enumerate(SPs, start = 1): #Extract the variable in the Stress Period
                sp = year[base:SP, :, :]
                base = SP
                #Sum the infiltration across the stress period
                if(units == 'inches'):
                    #Keep in inches (make the transformation later)
                    # Useful because otherwise ArcGIS can't read the values (too small)
                    sp = np.sum(sp, axis = 0) #inches
                elif(units == 'ms'):
                    #Transform into m/s
                    sp = np.sum(sp, axis = 0)*0.0254/(60*60*sp.shape[0]) #m/s
                else:
                    return print('Unrecognised unit. The available units are:\
                                 inches, ms (for meters/second)')
            
                #Save as Arc GRID ASCII file
                if(outpath != 'none'):
                    name = variable if name == 'name' else name
                    fname = f'{outpath}/{name}_{y}_SP{i}_inch.asc'
                    sp = pd.DataFrame(sp)
                    save_ArcGRID(sp, fname, xll, yll, size, self.metadata['nodata_value'])
                #Save in the 3D variable
                var3d[k, :, :] = sp
                k += 1
            s = e #It will get the subsequent day
        # self.SPsum = var3d
        # return self.SPsum
        print('End of the procedure')
        if outpath != 'none': print(f'The ASCII files are saved in {outpath}')
        return var3d
    
    def close(self):
        #Close the netCDF file
        self.netCDF.close()