# -*- coding: utf-8 -*-
"""
EOBSobject class definition

@author: paolo
"""

from datetime import date, timedelta
import glob
import netCDF4 as nc
import numpy as np
import os
import pandas as pd

class EOBSobject():
    
    def __init__(self, var, inpath = None, outpath = None, outname = None,
                 fname = None, API = False, folder = False, swb2 = False):
        """
        
        #var: name of the variable inside the netcdf file
        #folder: default is False (path to a single file). True, path to a folder
        
        
        Parameters
        ----------
        var : TYPE
            DESCRIPTION.
        inpath : TYPE, optional
            DESCRIPTION. The default is None.
        outpath : TYPE, optional
            DESCRIPTION. The default is None.
        outname : TYPE, optional
            DESCRIPTION. The default is None.
        fname : TYPE, optional
            DESCRIPTION. The default is None.
        API : TYPE, optional
            DESCRIPTION. The default is False.
        folder : TYPE, optional
            DESCRIPTION. The default is False.
        swb2 : TYPE, optional
            DESCRIPTION. The default is False.
        """
        #Store the info
        self.info = {
            'var': var,
            'for_swb2': swb2,
            'API': API
            }
        self.set_outname(outname)
        self.set_fname(fname)
        #Store the input path
        if API:
            pass
        elif folder: self.paths = { 'inpath': self.find_path(inpath, var) }
        else: self.paths = { 'inpath': inpath }
        self.set_outpath(outpath, folder)
    
    def load(self):
        if not self.info['API']:
            self.netcdf = nc.Dataset(self.paths['inpath'])
            #Store the units
            self.info['units'] = self.netcdf[self.info['var']].units
            self.info['missing_value'] = self.netcdf[self.info['var']]._FillValue
        else:
            #Run the API and load the desired dataset
            #May need additional parameters, as the version of E-OBS and such
            pass
    
    def print_metadata(self):
        #Raw metadata
        print('These are the original E-OBS metadata:')
        print(self.netcdf)
        #Print custom metadata
        # print('These are selected E-OBS metadata:')
        
    #---------------------------------------------------------
    #Perform operations on the original content
    
    def cut_space(self, coord, autosave = True, ext = False,
                 loncol = 'lon', latcol = 'lat', contourcell = 0,
                 saveformat = 'netcdf', readme = False):
        """
        coord: extremes of desired area
            provided as a pandas dataframe with loncol as the column containing
            longitude and latcol as the column containing latitude
        contourcell: number of contour cells to extract around the provided coordinates
        """
        la = self.netcdf['latitude'][:]
        lo = self.netcdf['longitude'][:]
        tool = round(la[1] - la[0], 1) * contourcell
        minlon = min(coord[loncol]) - tool
        minlat = min(coord[latcol]) - tool
        maxlon = max(coord[loncol]) + tool
        maxlat = max(coord[latcol]) + tool
        idx_lat = np.intersect1d(np.where(la > minlat), np.where(la < maxlat))
        idx_lon = np.intersect1d(np.where(lo > minlon), np.where(lo < maxlon))
        res = {
            'start_day': self.get_dates()[1],
            'end': self.get_dates()[2],
            'time': None,
            'idx_time': 0,
            'idx_lat': idx_lat,
            'idx_lon': idx_lon,
            'option': None
            }
        if autosave:
            if saveformat == 'netcdf': self.save_netcdf(res, method = 'cut_space', readme = readme)
            elif saveformat == 'arcgrid': self.save_arcgrid(res, method = 'cut_space', readme = readme)
        if ext: return res
    
    def cut_time(self, start, end, autosave = True, ext = False,
                 option = 'singleyear', day = False, saveformat = 'netcdf',
                 readme = False):
        """
        start, end: years (int) if day = False
           if day = True, they have to be in datetime.date format, ex: date(2014, 7, 20)
           day = True works only for option = 'bundle'
        option:
         - 'singleyear': single files, one for each year
         - 'bundle': one single file between the selected dates
        """
        #Get the E-OBS time limits from its metadata
        o, sd, ed = self.get_dates()
        #Create an array of the data real-world dates
        dt = pd.date_range(start = sd, end = ed)
        yR = dt.year
        yU = yR.unique()
        #Create an array of real-world dates from the origin of the dataset
        t = pd.date_range(start = o, end = ed)
        tyR = t.year
        if option == 'singleyear':
            for year in yU[np.where((yU >= start) & (yU <= end))]:
                idx_time = np.where(yR == year)[0]
                time = np.where(tyR == year)[0]
                if self.info['for_swb2']: time = self.transf_eobstime(time)
                res = {
                    'start_day': f"01/01/{year}",
                    'end': end,
                    'time': time,
                    'idx_time': idx_time,
                    'idx_lat': 0,
                    'idx_lon': 0,
                    'option': option
                    }
                if autosave: self.save_netcdf(res, method = 'cut_time', readme = readme)
                if ext: return res
        elif option == 'bundle':
            if day:
                idx_time = [*range(np.where(dt.date == start)[0].item() - 1, np.where(dt.date == end)[0].item() + 1)]
                idx_time = np.array(idx_time, dtype = np.int64)
                time = [*range(np.where(t.date == start)[0].item() - 1, np.where(t.date == end)[0].item() + 1)]
                time = np.array(time, dtype = np.int64)
            else:
                idx_time = np.where((yR >= start) & (yR <= end))[0]
                time = np.where((tyR >= start) & (tyR <= end))[0]
            res = {
                'start_day': f"{start.day}/{start.month}/{start.year}" if day else f"01/01/{start}",
                'end': end,
                'time': time,
                'idx_time': idx_time,
                'idx_lat': 0,
                'idx_lon': 0,
                'option': option
                }
            if autosave:
                if saveformat == 'netcdf': self.save_netcdf(res, method = 'cut_time', readme = readme)
                elif saveformat == 'arcgrid': self.save_arcgrid(res, method = 'cut_time', readme = readme)
            if ext: return res
        else:
            print('Wrong option inserted')
            return
    
    def cut_spacetime(self, coord, start, end, autosave = True, ext = False,
                      loncol = 'lon', latcol = 'lat', contourcell = 0,
                      option = 'singleyear', day = False, saveformat = 'netcdf',
                      readme = False):
        res_cs = self.cut_space(coord, False, True, loncol, latcol, contourcell)
        if option == 'singleyear':
            for year in range(start, end + 1):
                res_ct = self.cut_time(year, year, autosave = False, ext = True)
                res = {
                    'start_day': res_ct['start_day'],
                    'end': end,
                    'time': res_ct['time'],
                    'idx_time': res_ct['idx_time'],
                    'idx_lat': res_cs['idx_lat'],
                    'idx_lon': res_cs['idx_lon'],
                    'option': option
                    }
                if autosave:
                    if saveformat == 'netcdf': self.save_netcdf(res, method = 'cut_spacetime', readme = readme)
                    elif saveformat == 'arcgrid': self.save_arcgrid(res, method = 'cut_spacetime', readme = readme)
                if ext: return res
        elif option == 'bundle':
            res_ct = self.cut_time(start, end, False, True, option, day)
            res = {
                'start_day': res_ct['start_day'],
                'end': end,
                'time': res_ct['time'],
                'idx_time': res_ct['idx_time'],
                'idx_lat': res_cs['idx_lat'],
                'idx_lon': res_cs['idx_lon'],
                'option': option
                }
            if autosave:
                if saveformat == 'netcdf': self.save_netcdf(res, method = 'cut_spacetime', readme = readme)
                elif saveformat == 'arcgrid': self.save_arcgrid(res, method = 'cut_spacetime', readme = readme)
            if ext: return res
        else:
            print('Wrong option inserted')
            return
    
    #---------------------------------------------------------
    #Export the results

    def save_netcdf(self, res = None, method = 'raw', readme = True):
        #Check the 	Climate and Forecast Metadata Convention v1.4 (CF-v1.4)
        #and try to keep the metadata as they are defined there
        if readme: self.write_readme(res, method, savedas = 'NetCDF')
        outname = self.info['outname']
        if method == 'raw':
            #Saves the same dataset, just by applying a custom format
            _, tool, _ = self.get_dates()
            start_day = f"{tool.day}/{tool.month}/{tool.year}"
            la = self.get_lat(method)
            lo = self.get_lon(method)
            tout = np.ma.getdata(self.netcdf['time'][:])
        else:
            df = self.get_var(method, res['idx_time'], res['idx_lat'], res['idx_lon'])
            la = self.get_lat(method, res['idx_lat'])
            lo = self.get_lon(method, res['idx_lon'])
        
        if method == 'cut_space':
            description = "Clip in space of the E-OBS dataset"
            _, tool, _ = self.get_dates()
            start_day = f"{tool.day}/{tool.month}/{tool.year}"
            tout = np.ma.getdata(self.netcdf['time'][:])
        elif method == 'cut_time':
            description = "Clip in time of the E-OBS dataset"
            start_day = res['start_day']
            tout = res['time']
        elif method == 'cut_spacetime':
            description = "Clip in space and time of the E-OBS dataset"
            start_day = res['start_day']
            tout = res['time']
            
        if self.info['for_swb2']:
            la[::-1].sort()
            df = np.flip(df, axis = 1)
            df = np.around(df, 1)
            tunits = 'days since 1980-01-01 00:00:00 UTC'
            
            #Experiment to uniform also the files for SWB2:
            # - try to keep the same time enumeration as E-OBS and run SWB2
            # - try to keep the "start_day" as YYYY-MM-DD
        else:
            tunits = self.netcdf['time'].units
        
        fname = self.write_fname(method, res)
        #Create the netcdf dataset: netcdf3 for swb2, netcdf4 in general
        if self.info['for_swb2']: ds = nc.Dataset(fname, 'w', format = "NETCDF3_CLASSIC")
        else: ds = nc.Dataset(fname, 'w', format = 'NETCDF4')
        
        ## General metadata
        ds.description = description
        version = self.netcdf.__dict__[self.get_keys(self.netcdf.__dict__)[0]]
        ds.source = f"E-OBS {version}"
        ds.start_day = start_day
        #Syntax of start_day can be changed to YYYY-MM-DD,
        # but first I need to try run SWB2 in that configuration
        # ds.author = "paolocolombo1996@gmail.com"
        ds.reference_system = "WGS84"
        ds.proj4_string = "+proj=lonlat +datum=WGS84 +no_defs"
        
        ## Dimensions
        ds.createDimension('x', len(lo))
        ds.createDimension('y', len(la))
        ds.createDimension('time', None)
        
        ## Variables
        #x(x)
        x = ds.createVariable('x', 'd', ('x'))
        x.units = 'degrees'
        x.long_name = 'x coordinate of projection'
        x.standard_name = 'projection_x_coordinate'
        #y(y)
        y = ds.createVariable('y', 'd', ('y'))
        y.units = 'degrees'
        y.long_name = 'y coordinate of projection'
        y.standard_name = 'projection_y_coordinate'
        #Time
        time = ds.createVariable('time', 'd', ('time'))
        time.units = tunits
        #Day of the year
        yearday = ds.createVariable('yearday', 'h', ('time'))
        #Value
        value = ds.createVariable(outname, 'f4', ('time','y','x'), fill_value = -9999)
        value.units = self.info['units']
        value.missing_value = self.info['missing_value']
        value.coordinates = 'lat lon'
        
        ## Fill the variables
        x[:] = lo
        y[:] = la
        time[:] = tout
        yearday[:] = range(1, 366)
        value[:] = df if method != 'raw' else np.ma.getdata(self.netcdf[self.info['var']][:, :, :])
        
        #Close the file
        ds.close()

    def save_arcgrid(self, res = None, method = 'raw', createfolder = True,
                     custom = False, customdf = None, customname = None,
                     readme = False):
        """
        Save the E-OBS dataset as daily ArcGRID files
        """
        if readme: self.write_readme(res, method, savedas = 'ArcGRID')
        if method == 'raw':
            la = self.get_lat(method)
            lo = self.get_lon(method)
            t = self.get_time(method)
        else:
            if not custom:
                df = self.get_var(method, res['idx_time'], res['idx_lat'], res['idx_lon'])
                df = np.flip(df, axis = 1)
                t = self.get_time(method, res)
            else:
                t = [0 for i in range(customdf.shape[0])]
            la = self.get_lat(method, res['idx_lat'])
            lo = self.get_lon(method, res['idx_lon'])
        
        namefolder = self.info['outname']
        namefile = namefolder.upper() if self.info['for_swb2'] else namefolder
        if (namefolder == 'precip') & (self.info['for_swb2']): namefile = 'PRCP'
        outpath = self.paths['outpath']
        if createfolder:
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            outpath = f"{outpath}/{namefolder}"
            if not os.path.exists(outpath):
                os.makedirs(outpath)
        
        size = round(la[1] - la[0], 1)
        xll = lo[0] - size/2
        yll = la[0] - size/2
        nodata = self.info['missing_value']
        
        for i in range(len(t)):
            if not custom:
                y, m, d = self.transf_eobsdate(t[i], number = True)
                fname = f'{outpath}/{namefile}_{y}_{m}_{d}.asc'
                if method != 'raw':
                    tool = round(pd.DataFrame(df[i, :, :]), 1)
                else:
                    tool = round(pd.DataFrame(np.flip(np.ma.getdata(self.netcdf[self.info['var']][i, :, :]), 0)), 1)
            else:
                tool = pd.DataFrame(np.flip(customdf[i, :, :], 0))
                fname = f'{outpath}/{customname}_{i+1}.asc'
            tool.to_csv(fname, sep = ' ', header = False, index = False)
            header = f'ncols         {tool.shape[1]}\nnrows         {tool.shape[0]}\nxllcorner     {xll}\nyllcorner     {yll}\ncellsize      {size}\nNODATA_value  {nodata}'
            with open(fname, 'r+') as f:
                content = f.read()
                f.seek(0, 0)
                f.write(header + '\n' + content)   

    #----------------------------------------------------------
    #Generation of statistics or additional information
    
    def SP_sum(self, SPs, checkleap = True, export = True, store = True,
               units = None, method = 'raw', coord = None, start = None, end = None,
               loncol = 'lon', latcol = 'lat', contourcell = 0, option = 'bundle'):
        '''
        Sum of the variable in provided "stress periods" (SPs)
        '''
        
        if method == 'cut_space':
            res = self.cut_space(coord, False, True, loncol, latcol, contourcell = 0)
        elif method == 'cut_time':
            res = self.cut_time(start, end, False, True, option)
        elif method == 'cut_spacetime':
            res = self.cut_spacetime(coord, start, end, False, True, 
                                     loncol, latcol, contourcell, option)
        if method == 'raw':
            df = np.ma.getdata(self.netcdf[self.info['var']][:, :, :])
        else:
            df = self.get_var(method, res['idx_time'], res['idx_lat'], res['idx_lon'])
        
        SPs = np.cumsum(SPs)
        units = units if units else self.info['units']
        if start is None: start = self.get_dates()[1].year
        if end is None: end = self.get_dates()[2].year
        period = range(start, end+1)
        s, e, k = 0, 0, 0
        #Create the 3D variable
        var3d = np.zeros((len(period)*len(SPs), df.shape[1], df.shape[2]))
        
        for y in period:
            #Extract a single year
            e += self.leap(y) if checkleap else 365
            year = df[s:e, :, :]
            #Set up a counter
            base = 0
            for i, SP in enumerate(SPs, start = 1):
                if (checkleap) & (self.leap(y) == 366): SP = SP+1 #& (i == 1)
                #Extract the variable in the Stress Period
                sp = year[base:SP, :, :]
                base = SP
                #Sum the variable in the stress period
                if units == self.info['units']:
                    #Keep the E-OBS original units (mm)
                    sp = np.sum(sp, axis = 0) #mm
                elif units == 'ms':
                    #Transform from mm into m/s
                    sp = np.sum(sp, axis = 0)*1000/(60*60*24*sp.shape[0]) #m/s
                else:
                    return print('Unrecognised unit. The available units are:\
                                 ms (for meters/second)')
                #Save in the 3D variable
                var3d[k, :, :] = sp
                k += 1
            s = e #It will get the subsequent day
        if export: self.save_arcgrid(res, method, False, True, var3d, f'E-OBS_SP_sum_{units}')
        if store: self.SP_sum_df = var3d
    
    def write_readme(self, res = None, method = None, savedas = None, path = None, outname = None):
        if not path: path = self.paths['outpath']
        if not outname: outname = self.info['outname']
        
        readme = open(f"{path}/readme_{outname}.txt", "w")
        header = f"Source of the data: {self.netcdf.ncattrs()[0]} {getattr(self.netcdf, self.netcdf.ncattrs()[0])}\n"
        readme.write(header)
        o, s, e = self.get_dates()
        info = f"Original dates: \n \
    Origin: {o}\n \
    First data available: {s}\n \
    Last data available: {e}\n"
        readme.write(info)
        if method:            
            if method == 'cut_space':
                op = 'Crop in space'
            elif method == 'cut_time':
                op = 'Crop in time'
            elif method == 'cut_spacetime':
                op = 'Crop in both time and space' 
            elif method == 'raw':
                op = 'No performed operations'
            readme.write(f"Performed operation: {op}\n")
        if savedas: readme.write(f"Output format: {savedas}")
        if res:
            nrow = len(res['idx_lat']) if (method) and (method != 'cut_time') else self.netcdf['latitude'].shape[0]
            ncol = len(res['idx_lon']) if (method) and (method != 'cut_time') else self.netcdf['longitude'].shape[0]
            readme.write(f"Resulting data information:\n \
    Starting date: {res['start_day']}\n \
    Ending date (end of year): {res['end']}\n \
    Number of rows: {nrow}\n \
    Number of columns: {ncol}\n")
        readme.close()
    
    #----------------------------------------------------------
    #General operations
    
    def close_netcdf(self):
        self.netcdf.close()
    
    def find_path(self, inpath, var, ext = '.nc'):
        fls = glob.glob(f'{inpath}/*{ext}')
        for fl in fls:
            varin = fl[len(inpath)+1:].split('_', 1)[0]
            if var == varin:
                pos =  fls.index(fl)
                break
        return fls[pos]
    
    def get_dates(self):
        """
        Returns
        -------
        origin : datetime.date
            Origin of the progressive number used as time inside E-OBS
        start : datetime.date
            Date of the first data available inside the E-OBS file
        end : datetime.date
            Date of the last data available inside the E-OBS file
        """
        origin = self.netcdf['time'].units.split(" ")[2]
        origin = date.fromisoformat(origin)
        ndays = np.ma.getdata(self.netcdf['time'][0]).item()
        start = origin + timedelta(ndays)
        ndays = np.ma.getdata(self.netcdf['time'][-1]).item()
        end = origin + timedelta(ndays)
        return origin, start, end
    
    def get_keys(self, dict):
        """
        Returns the dictionary keys of a dictonary as a list
        Got from:
        https://www.geeksforgeeks.org/python-get-dictionary-keys-as-a-list/
        """
        #Other method
        # return list(dict.keys())
        return [*dict]
    
    def get_lat(self, method, idx_lat = None):
        if (method == 'cut_space') or (method == 'cut_spacetime'): 
            return self.netcdf['latitude'][idx_lat]
        else:
            return self.netcdf['latitude'][:]
    
    def get_lon(self, method, idx_lon = None):
        if (method == 'cut_space') or (method == 'cut_spacetime'): 
            return self.netcdf['longitude'][idx_lon]
        else:
            return self.netcdf['longitude'][:]
    
    def get_time(self, method, res = None):
        if (method == 'cut_time') or (method == 'cut_spacetime'):
            return res['time']
        else:
            return self.netcdf['time'][:]
    
    def get_var(self, method, idx_time = None, idx_lat = None, idx_lon = None):
        if method == 'raw':
            return np.ma.getdata(self.netcdf[self.info['var']][:, :, :])
        elif method == 'cut_space':
            return np.ma.getdata(self.netcdf[self.info['var']][:, idx_lat, idx_lon])
        elif method == 'cut_time':
            return np.ma.getdata(self.netcdf[self.info['var']][idx_time, :, :])
        elif method == 'cut_spacetime':
            return np.ma.getdata(self.netcdf[self.info['var']][idx_time, idx_lat, idx_lon])
    
    def leap(self, y):
        #input: year (int)
        #output: number of days (int)
        if((y%4 == 0) | (y%400 == 0)):
            return 366
        else:
            return 365
    
    def set_fname(self, fname):
        self.info['fname'] = fname
    
    def set_outname(self, outname):
        self.info['outname'] = outname if outname else self.info['var']
    
    def set_outpath(self, outpath, folder):
        tool = self.paths['inpath']
        if not folder:
            tool = tool.split('/')[:-1]
            tool = '/'.join(tool)
        else: tool = tool.split('\\')[0]
        self.paths['outpath'] = outpath if outpath else tool
    
    def set_timeunit(self, units = None):
        if not units:
            units = self.netcdf['time']['units'] #check if this works
        self.units['time'] = units
    
    def transf_eobstime(self, x, to = date(1980, 1, 1)):
        """
        Returns the time as number of days starting from a desired point in time
        
        Parameters
        ----------
        x : int
            date in the E-OBS format
        to : datetime.date
            time reference to be transformed to. Default is Daymet,
            starting from 1980-01-01
        """
        dstart = to
        estart = date(1950, 1, 1)
        k = dstart - estart
        #For SWB2, add 0.5
        y = x - k.days + 0.5 if self.info['for_swb2'] else x - k.days
        return y
    
    def transf_eobsdate(self, x, number = False):
        """
        Returns the year, month and day corresponding to the number given (x)
        If x is a plain number (not a variable), "number" must be set to True
        """
        from datetime import date, timedelta
        start = date(1980, 1, 1) if self.info['for_swb2'] else date(1950, 1, 1)
        days = x.item() if number else x
        end = start + timedelta(days = days)
        return end.year, end.strftime('%m'), end.strftime('%d')
    
    def write_fname(self, method, res):
        """
        Writes fname for save_netcdf function, starting from the method utilized        
        """
        if self.info['fname']: return self.info['fname']
        outpath = self.paths['outpath']
        outname = self.info['outname']
        fname = f'{outpath}/{outname}_EOBS_{method}'
        if (method == 'cut_time' or method == 'cut_spacetime') and res['option'] != 'bundle':
            tool = res['start_day'].split('/')[2]
            fname = f'{fname}_{tool}'
        elif (method == 'cut_time' or method == 'cut_spacetime') and res['option'] == 'bundle':
            tool = f"{res['start_day'].split('/')[2]}_{res['end']}"
            fname = f'{fname}_{tool}'
        fname = f'{fname}.nc'
        return fname