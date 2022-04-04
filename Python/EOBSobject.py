# -*- coding: utf-8 -*-
"""
EOBSobject class definition

@author: paolo
"""

import glob
import netCDF4 as nc
import numpy as np
import pandas as pd
from datetime import date, timedelta

class EOBSobject():
    
    def __init__(self, inpath, var, outpath = 'none', outname = 'none',
                 folder = True, swb2 = False):
        #Store the info
        self.info = {
            'var': var,
            'for_swb2': swb2
            }
        self.set_outname(outname)
        #Store the path
        #folder condition: default is True, path to a folder
        #if path to a single file
        if folder: self.paths = { 'inpath': self.find_path(inpath, var) }
        else: self.paths = { 'inpath': inpath }
        self.set_outpath(outpath, folder)
        #Store the units
        self.units = {
            # 'variable': #save here the unit of the variable
            }
    
    def load(self):
        self.netcdf = nc.Dataset(self.paths['inpath'])
    
    def print_metadata(self):
        #Raw metadata
        print('These are the original E-OBS metadata:')
        print(self.netcdf)
        #Print custom metadata
        # print('These are selected E-OBS metadata:')
        
    #---------------------------------------------------------
    #NETCDF section
    
    def cut_space(self, coord, save = True, internal = False,
                 loncol = 'lon', latcol = 'lat', contourcell = 0):
        #coord: extremes of desired area
        #contourcell: number of contour cells to extract around the provided
        # coordinates
        
        la = self.netcdf['latitude'][:]
        lo = self.netcdf['longitude'][:]
        tool = round(la[0] - la[1], 1) * contourcell
        minlon = min(coord[loncol]) - tool
        minlat = min(coord[latcol]) - tool
        maxlon = max(coord[loncol]) + tool
        maxlat = max(coord[latcol]) + tool
        idx_lat = np.intersect1d(np.where(la > minlat), np.where(la < maxlat))
        idx_lon = np.intersect1d(np.where(lo > minlon), np.where(lo < maxlon))
        res = {
            'idx_time': 0,
            'idx_lat': idx_lat,
            'idx_lon': idx_lon
            }
        if internal: return res
        if save: self.save_netcdf(res, 'cut_space')
    
    def cut_time(self, start, end, save = True, internal = False,
                 option = 'singleyear', day = False):
        #start, end: years (int) if day = False
        #   if day = True, they have to be in datetime.date format, ex: date(2014, 7, 20)
        #   day = True works only for option = 'bundle'
        #option:
        # - 'singleyear': single files, one for each year
        # - 'bundle': one single file between the selected dates
        
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
                    'time': time,
                    'idx_time': idx_time,
                    'idx_lat': 0,
                    'idx_lon': 0,
                    'option': option
                    }
                if save: self.save_netcdf(res, method = 'cut_time')
                if internal: return res
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
                'time': time,
                'idx_time': idx_time,
                'idx_lat': 0,
                'idx_lon': 0,
                'option': option
                }
            if save: self.save_netcdf(res, method = 'cut_time')
            if internal: return res
        else:
            print('Wrong option inserted')
            return
    
    def cut_spacetime(self, coord, start, end, save = True, internal = False,
                      loncol = 'lon', latcol = 'lat', contourcell = 0,
                      option = 'singleyear', day = False):
        res_ca = self.cut_space(coord, False, True,
                                loncol, latcol, contourcell)
        if option == 'singleyear':
            for year in range(start, end + 1):
                res_ct = self.cut_time(year, year, save = False, internal = True)
                res = {
                    'start_day': res_ct['start_day'],
                    'time': res_ct['time'],
                    'idx_time': res_ct['idx_time'],
                    'idx_lat': res_ca['idx_lat'],
                    'idx_lon': res_ca['idx_lon'],
                    'option': option
                    }
                if save: self.save_netcdf(res, method = 'cut_spacetime')
                if internal: return res
        elif option == 'bundle':
            res_ct = self.cut_time(start, end, False, True, option, day)
            res = {
                'start_day': res_ct['start_day'],
                'time': res_ct['time'],
                'idx_time': res_ct['idx_time'],
                'idx_lat': res_ca['idx_lat'],
                'idx_lon': res_ca['idx_lon'],
                'option': option
                }
            if save: self.save_netcdf(res, method = 'cut_spacetime')
            if internal: return res
        else:
            print('Wrong option inserted')
            return
    
    def save_netcdf(self, res = None, method = 'raw'):
        # outpath = self.paths['outpath']
        outname = self.info['outname']
        # check = False
        
     
        
        #ottimizzare il passaggio di variabili tra funzioni
        #df può anche non essere passato, passando solo i vari idx
        #e estraendolo direttamente qui
        
        if method == 'raw':
            #Saves the same dataset, just by applying a custom format
            _, tool, _ = self.get_dates()
            start_day = f"{tool.day}/{tool.month}/{tool.year}"
            df = self.get_var(method)
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
            description = "Clip in sapce and time of the E-OBS dataset"
            start_day = res['start_day']
            tout = res['time']
            
        if self.info['for_swb2']:
            la[::-1].sort()
            df = np.flip(df, axis = 1)
            df = np.around(df, 1)
            tunits = 'days since 1980-01-01 00:00:00 UTC'
            #Expertiment:
            #   try to keep the same numeration as E-OBS
        else:
            tunits = self.netcdf['time'].units
        
        fname = self.write_fname(method, res)
        
        # if not check: fname = f'{outpath}/{outname}_EOBS_{method}.nc'
        if self.info['for_swb2']: ds = nc.Dataset(fname, 'w', format = "NETCDF3_CLASSIC")
        else: ds = nc.Dataset(fname, 'w')
        
        ## General metadata
        ds.description = description
        ds.source = "E-OBS v24.0"
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
        value.missing_value = -9999.0
        value.coordinates = 'lat lon'
        
        ## Fill the variables
        x[:] = lo
        y[:] = la
        time[:] = tout
        yearday[:] = range(1, 366)
        value[:] = df
        
        #Close the file
        ds.close()
    
    def save_arcgrid(self, code = 0):
        self.paths['outpath']
        if code == 0: self.cut_space() #call cut_space
        #call cut_time
        #call cut_spacetime
    
    #----------------------------------------------------------
    #ASCII section
    
    #anzi, basta definire bene le funzioni cut in modo che ritornino qualcosa
    #utilizzabile sia da una funzione save_netcdf sia da una save_arcgrid
    
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
        origin = self.netcdf['time'].units.split(" ")[2]
        origin = date.fromisoformat(origin)
        ndays = np.ma.getdata(self.netcdf['time'][0]).item()
        start = origin + timedelta(ndays)
        ndays = np.ma.getdata(self.netcdf['time'][-1]).item()
        end = origin + timedelta(ndays)
        return origin, start, end
    
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
    
    def get_var(self, method, idx_time = None, idx_lat = None, idx_lon = None):
        if method == 'raw':
            df = np.ma.getdata(self.netcdf[self.info['var']][:, :, :])
        elif method == 'cut_time':
            df = np.ma.getdata(self.netcdf[self.info['var']][idx_time, :, :])
        elif method == 'cut_space':
            df = np.ma.getdata(self.netcdf[self.info['var']][:, idx_lat, idx_lon])
        elif method == 'cut_spacetime':
            df = np.ma.getdata(self.netcdf[self.info['var']][idx_time, idx_lat, idx_lon])
        return df
    
    def set_outname(self, outname):
        self.info['outname'] = outname if outname != 'none' else self.info['var']
    
    def set_outpath(self, outpath, folder):
        tool = self.paths['inpath']
        if not folder:
            tool = tool.split('/')[:-1]
            tool = '/'.join(tool)
        else: tool = tool.split('\\')[0]
        self.paths['outpath'] = outpath if outpath != 'none' else tool
    
    def set_timeunit(self, units = 'none'):
        if units == 'none':
            units = self.netcdf['time']['units'] #check if this works
        self.units['time'] = units
    
    def transf_eobstime(self, y, to = date(1980, 1, 1)):
        #to: time reference to be transformed to
        #   standard: Daymet, starting from 1980-01-01
        dstart = to
        estart = date(1950, 1, 1)
        k = dstart - estart
        #For SWB2, add 0.5
        y1 = y - k.days + 0.5 if self.info['for_swb2'] else y - k.days
        return y1
    
    def write_fname(self, method, res):
        outpath = self.paths['outpath']
        outname = self.info['outname']
        fname = f'{outpath}/{outname}_EOBS_{method}'
        if (method == 'cut_time' or method == 'cut_spacetime') and res['option'] != 'bundle':
            tool = res['start_day'].split('/')[2]
            fname = f'{fname}_{tool}'
        fname = f'{fname}.nc'
        return fname
        
        
        
        