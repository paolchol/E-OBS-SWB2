# -*- coding: utf-8 -*-
"""
Class "SWB2output" for SWB2 netCDF output handling
It can be called using:
from Python.SWB2output import SWB2output
The working directory has to be set in ./E-OBS-SWB2 for this to work

@author: paolo
"""

import netCDF4 as nc
import numpy as np
import pandas as pd

class SWB2output():
    
    def __init__(self, path):
        """
        Load the NetCDF provided by SWB2 as output and extract its main 
        metadata

        Parameters
        ----------
        path : str
            Path to the SWB2 NetCDF output file.
        """        
        self.path = path
        self.netCDF = nc.Dataset(path)
        start, end = self.getdates(self.netCDF.variables['time'].units,
                              self.netCDF.variables['time'].shape[0])
        self.metadata = {
            "start_date": start,
            "end_date": end,
            # "variables": self.getkeys(self.netCDF.variables),
            # "main_variable": self.getkeys(self.netCDF.variables)[3], 
            # "units": self.netCDF[self.getkeys(self.netCDF.variables)[3]].units,
            # "nodata_value": self.netCDF[self.getkeys(self.netCDF.variables)[3]]._FillValue
            "variables": [*self.netCDF.variables],
            "main_variable": [*self.netCDF.variables][3], 
            "units": self.netCDF[[*self.netCDF.variables][3]].units,
            "nodata_value": self.netCDF[[*self.netCDF.variables][3]]._FillValue
            }
        self.results = {}
    
    def print_metadata(self):
        """
        Print original and custom metadata
        """
        print('Original metadata')
        print(self.netCDF)
        print('Extracted metadata')
        print(self.metadata)
    
    #Operations on the data
    #----------------------
    
    def extract(self):
        """
        Returns
        -------
        numpy.array
            Numpy 3D array containing the main variable of SWB2 output provided.
        """
        variable = self.metadata['main_variable']
        return np.ma.getdata(self.netCDF[variable][:, :, :])

    def SP_sum_total(self, SPs, units='ms', retval=False):
        """
        Sum SWB2 net infiltration over stress periods defined for the
        entire simulation.

        SPs: list of stress-period lengths in days.
        """

        variable = 'net_infiltration'

        print(f'Performing the sum of {variable} over the stress periods')
        print(f'Output unit measure: {units}')

        data = np.ma.getdata(self.netCDF[variable])

        print("Variable:", variable)
        print("Data shape:", data.shape)

        if data.ndim != 3:
            raise ValueError(
                f"Expected SWB2 variable to have 3 dimensions "
                f"(time, y, x), but got shape {data.shape}"
            )

        n_days = data.shape[0]

        print("Number of daily records:", n_days)
        print("Number of SPs:", len(SPs))
        print("SP lengths:", SPs)
        print("Total SP days:", sum(SPs))

        if sum(SPs) != n_days:
            raise ValueError(
                f"Stress periods cover {sum(SPs)} days, "
                f"but SWB2 output contains {n_days} days."
            )

        # Create output
        var3d = np.zeros(
            (len(SPs), data.shape[1], data.shape[2]),
            dtype=np.float32
        )

        start = 0

        for i, SP in enumerate(SPs):

            end = start + SP

            # Extract stress period
            sp = data[start:end, :, :]

            # Sum daily recharge
            if units == 'inches':
                # Total recharge in inches
                result = np.sum(sp, axis=0)

            elif units == 'ms':
                # inches/day -> m/s
                result = (
                    np.sum(sp, axis=0)
                    * 0.0254
                    / (86400 * SP)
                )

            else:
                raise ValueError(
                    "Unrecognised unit. Available units: "
                    "'inches', 'ms'"
                )

            var3d[i, :, :] = result

            print(
                f"SP{i+1}: days {start}–{end-1} "
                f"({SP} days)"
            )

            start = end

        self.results['SPsum3d'] = var3d

        if retval:
            return var3d
        
    def sumtot(self, outpath = None, name = 'name'):
        """
        Returns the sum of the main variable over the whole time period

        Parameters
        ----------
        outpath : str, optional
            Path fro the output file. The default is None.
        name : str, optional
            First part of the output file name. The default is 'name'.

        Returns
        -------
        pandas.DataFrame
            Dataframe of the sum performed over the whole time period.
        """
        variable = self.metadata["main_variable"]
        self.sumtotdf = np.sum(np.ma.getdata(self.netCDF[variable][:,:,:]), axis = 0)*0.0254 #meters
        
        if outpath:
            df = pd.DataFrame(self.sumtotdf)
            size = round(np.ma.getdata(self.netCDF['x'][1]).item() - np.ma.getdata(self.netCDF['x'][0]).item())
            xll = round(np.ma.getdata(self.netCDF['x'][0]).item()) - size/2
            yll = round(np.ma.getdata(self.netCDF['y'][-1]).item()) - size/2
            fname = f'{outpath}/{name}_{variable}_sum_tot.asc'
            self.save_ArcGRID(df, fname, xll, yll, size, -9999)
            print(f'ArcGRID saved in {outpath} as: {name}_{variable}_sum_tot.asc')
        return self.sumtotdf
    
    def SP_sum(self, SPs, outpath = 'none', name = 'name', units = 'none',
               retval = False, checkleap = True):
        """
        Perform a sum of the main variable over the stress periods provided
        as SPs
        
        Returns a 3D variable containing all the sums as the 0 index
        If outpath is specified, saves an ArcASCII GRID file for each sum
        """
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
        
        for y in period:
            #Extract a single year
            e += self.leap(y) if checkleap else 365
            year = np.ma.getdata(self.netCDF[variable][s:e, :, :])
            #Set up a counter
            base = 0
            for i, SP in enumerate(SPs, start = 1):
                if (checkleap) & (self.leap(y) == 366): SP = SP+1 #& (i == 1)
                #Extract the variable in the Stress Period
                sp = year[base:SP, :, :]
                base = SP
                #Sum the variable in the stress period
                if(units == 'inches'):
                    #Keep in inches (make the transformation later)
                    #Useful because otherwise ArcGIS can't read the values (too small)
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
                    fname = f'{outpath}/{name}_{y}_SP{i}_{units}.asc'
                    sp = pd.DataFrame(sp)
                    self.save_ArcGRID(sp, fname, xll, yll, size, self.metadata['nodata_value'])
                #Save in the 3D variable
                var3d[k, :, :] = sp
                k += 1
            s = e #It will get the subsequent day
        
        print('End of the procedure')
        if outpath != 'none': print(f'The ASCII files are saved in {outpath}')
        self.results['SPsum3d'] = var3d
        if retval: return var3d
    
    def obtain_df_SP(self, index = 'index'):
        """
        Obtain a dataframe containing the results of the sum operated in each 
        stress period. Needs method SP_sum() to be performed before.
        The result is saved in self.results['dfSP'].
        
        Parameters
        ----------
        index : str, optional
            Name to be given to the index column in the resulting dataframe.
            The default is 'index'.
        """
        if 'SPsum3d' in self.results:
            df = pd.DataFrame(self.results['SPsum3d'][0, :, :])
        else:
            print('Error: Run SP_sum() method before running obtain_df_SP()')
            return
        df.insert(0, 'nrow', df.index.values)
        df = pd.melt(df, id_vars = 'nrow', var_name = 'ncol', value_name = 'SP1')
        df['nrow'] = df['nrow'] + 1
        df['ncol'] = df['ncol'] + 1
        df = self.insertind(df, df['nrow'], df['ncol'], name = index)
        
        for i in range(1, self.results['SPsum3d'].shape[0]):
            df = pd.DataFrame(self.results['SPsum3d'][i, :, :])
            df.insert(0, 'nrow', df.index.values)
            df = pd.melt(df, id_vars = 'nrow', var_name = 'ncol', value_name = f'SP{i+1}')
            #insert indicator
            if f'SP{i+1}' not in df.columns:
                df.insert(len(df.columns), f'SP{i+1}', df[f'SP{i+1}'])
        self.results['dfSP'] = df
    
    #Export functions
    #----------------
    
    def save_ArcGRID(self, df, fname, xll, yll, size, nodata):
        """
        Saves an ArcGRID from the provided dataframe

        Parameters
        ----------
        df : pandas.DataFrame
            The dataframe you want to save as an ArcGRID.
        fname : str
            The name of the output file.
        xll : float
            x coordinate of the left bottom corner (lon).
        yll : float
            y coordinate of the left bottom corner (lat).
        size : float
            Cell size (m)
        nodata : float
            Value assigned to nodata.
        """
        def line_prepender(filename, line):
            with open(filename, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                f.write(line + '\n' + content)
                #line.rstrip('\r\n') if you want ot remove something from line
        df.to_csv(fname, sep = ' ', header = False, index = False)
        header = f'ncols         {len(df.columns)}\nnrows         {len(df.index)}\nxllcorner     {xll}\nyllcorner     {yll}\ncellsize      {size}\nNODATA_value  {nodata}'
        line_prepender(fname, header)
    
    #General functions
    #-----------------   
    
    def close(self):
        """
        Close the netCDF file
        """
        self.netCDF.close()
    
    def getdates(self, string, n):
        """
        Get the starting and ending dates of an SWB2 netCDF output file
        """
        from datetime import date, timedelta
        datestr = string.split(' ')[2]
        y, m, d = datestr.split('-')
        start = date(int(y), int(m), int(d))
        end = start + timedelta(days = n-1)
        return start, end
    
    def getkeys(dict):
        """
        Returns a dictionary's keys as a list
        Got from:
        https://www.geeksforgeeks.org/python-get-dictionary-keys-as-a-list/
        Other method:
        return list(dict.keys())

        Parameters
        ----------
        dict : dict
            Dictionary from which extract the keys.
        
        Returns
        -------
        list
            Dictionary keys.
        """
        return [*dict]
    
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
    
    def leap(self, y):
        """
        Returns the number of days of the year provided

        Parameters
        ----------
        y : int
            Year for which to know the number of days.

        Returns
        -------
        int
            number of days of the year provided.
        """
        if (y%4 == 0) | (y%400 == 0):
            return 366
        else:
            return 365