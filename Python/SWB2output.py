# -*- coding: utf-8 -*-
"""
Class "SWB2output" for SWB2 netCDF output handling
It can be called using:
from Python.SWB2output import SWB2output

The working directory has to be set in ./E-OBS-SWB2 for this to work

@author: paolo
"""

# import os
import netCDF4 as nc
import numpy as np
import pandas as pd

#The directory has to be set in ./E-OBS-SWB2 for this to work
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
        self.results = {}
    
    def print_md(self):
        #Prints the original netCDF metadata
        print(self.netCDF)
    
    def extract(self):
        #Returns the main variable of the output as a numpy 3d array
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
            fname = f'{outpath}/{name}_{variable}_sum_tot.asc'
            save_ArcGRID(df, fname, xll, yll, size, -9999)
            print(f'ArcGRID saved in {outpath} as: {name}_{variable}_sum_tot.asc')
        return self.sumtotdf
    
    def SP_sum(self, SPs, outpath = 'none', name = 'name', units = 'none',
               retval = False):
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
                    sp = np.sum(sp, axis = 0)*0.0254/(60*60*24*sp.shape[0]) #m/s
                else:
                    return print('Unrecognised unit. The available units are:\
                                 inches, ms (for meters/second)')
            
                if(outpath != 'none'):
                    #Save as Arc GRID ASCII file
                    name = variable if name == 'name' else name
                    fname = f'{outpath}/{name}_{y}_SP{i}_inch.asc'
                    sp = pd.DataFrame(sp)
                    save_ArcGRID(sp, fname, xll, yll, size, self.metadata['nodata_value'])
                #Save in the 3D variable
                var3d[k, :, :] = sp
                k += 1
            s = e #It will get the subsequent day
        
        print('End of the procedure')
        if outpath != 'none': print(f'The ASCII files are saved in {outpath}')
        self.results['SPsum3d'] = var3d
        if retval: return var3d
    
    def obtain_df_SP(self, index = 'none'):
        if 'SPsum3d' in self.results:
            df = pd.DataFrame(self.results['SPsum3d'][0, :, :])
        else:
            print('Error: Run SP_sum() method before running obtain_df_SP()')
            return
        df.insert(0, 'nrow', df.index.values)
        df = pd.melt(df, id_vars = 'nrow', var_name = 'ncol',
         value_name = 'SP1')
        df['nrow'] = df['nrow'] + 1
        df['ncol'] = df['ncol'] + 1
        df = self.insertind(df, df['nrow'], df['ncol'], name = index)
        
        for i in range(1, self.results['SPsum3d'].shape[0]):
            df = pd.DataFrame(self.results['SPsum3d'][i, :, :])
            df.insert(0, 'nrow', df.index.values)
            df = pd.melt(df, id_vars = 'nrow', var_name = 'ncol',
             value_name = f'SP{i+1}')
            if f'SP{i+1}' not in df.columns:
                df.insert(len(df.columns), f'SP{i+1}', df[f'SP{i+1}'])
        self.results['dfSP'] = df
    
    def insertind(self, df, r, c, pos = 0, name = 'none'):
        # name = self.info['id'] if name == 'none' else name
        newc = []
        for i in range(len(r)):
            newc += [f'{r[i]}X{c[i]}']
        if (name not in df.columns):
            df.insert(pos, name, newc)
        else:
            df[name] = newc
        return df
    
    def close(self):
        #Close the netCDF file
        self.netCDF.close()