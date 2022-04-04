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
                 folder = True, swb2 = True):
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
        self.netcdf
        #Print custom metadata
        
    #---------------------------------------------------------
    #NETCDF section
    
    def cut_area(self, coord, save = True, internal = False,
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
        # df = self.netcdf[self.info['var']][:, idx_lat, idx_lon]
        # df = np.ma.getdata(df)
        res = {
            # 'df': df,
            # 'la': la[idx_lat],
            # 'lo': lo[idx_lon],
            'idx_time': 0,
            'idx_lat': idx_lat,
            'idx_lon': idx_lon
            }
        if internal: return res
        if save: self.save_netcdf(res, 'cut_area')
    
    def cut_time(self, start, end, save = True, internal = False,
                 option = 'singleyear'):
        #option:
        # - 'singleyear': single files, one for each year
        # - 'complete': one single file between the selected dates
        
        #Get the E-OBS time limits from its metadata
        o, sd, ed = self.get_dates()
        dt = pd.date_range(start = sd, end = ed)
        yR = dt.year
        yU = yR.unique()
        #E-OBS date definition (to insert it in the output .nc file)
        t = pd.date_range(start = o, end = ed)
        tyR = t.year
        
        if option == 'singleyear':
            for year in yU[np.where((yU >= start) & (yU <= end))]:
                idx_time = np.where(yR == year)[0]
                # df = self.netcdf[self.info['var']][idx_time, : , :]
                time = np.where(tyR == year)[0]
                if self.info['for_swb2']: time = self.transf_daymettime(time)
                res = {
                    # 'df': df,
                    # 'la': self.netcdf['latitude'][:],
                    # 'lo': self.netcdf['longitude'][:],
                    'year': year,
                    'time': time,
                    'idx_time': idx_time,
                    'idx_lat': 0,
                    'idx_lon': 0
                    }
                if save: self.save_netcdf(res, method = 'cut_time')
                if internal: return res
                
                #serve for loop in cut_areatime
                #start e end uguali, incrementano di 1
                #oppure tagliare il tempo direttamente in cut_areatime
        # elif option == 'complete':
            
        else:
            print('Wrong option inserted')
            return        
    
    def cut_areatime(self, save = False):
        self.cut_area(save = False, internal = True)
        
        #returns res: take idx_lat and idx_len
        
        self.cut_time(save = False, internal = True)
        
        #returns idx_time
        
        if save: self.save_netcdf()
    
    def save_netcdf(self, res = None, method = 'raw'):
        outpath = self.paths['outpath']
        outname = self.info['outname']
        
     
        
        #ottimizzare il passaggio di variabili tra funzioni
        #df può anche non essere passato, passando solo i vari idx
        #e estraendolo direttamente qui
        

        
        if method == 'raw':
            #save the same dataset only in a new format
            
            df = self.get_var(method)
            la = self.get_lat(method)
            lo = self.get_lon(method)
            tout = np.ma.getdata(self.netcdf['time'][:])
        else:
            df = self.get_var(method, res['idx_time'], res['idx_lat'], res['idx_lon'])
            la = self.get_lat(method, res['idx_lat'])
            lo = self.get_lon(method, res['idx_lon'])
        
        if method == 'cut_area':
            description = "Clip in space of the E-OBS dataset"
            tout = np.ma.getdata(self.netcdf['time'][:])
        elif method == 'cut_time':
            description = "Clip in time of the E-OBS dataset"
            start_day = f"01/01/{res['year']}"
            # end_day = calcola l'end day
            tout = res['time']
        elif method == 'cut_areatime':
            tout = res['time']
            pass
            
        if self.info['for_swb2']:
            la[::-1].sort()
            df = np.flip(df, axis = 1)
            df = np.around(df, 1)
            tunits = 'days since 1980-01-01 00:00:00 UTC' #try to keep the same numeration as E-OBS
        else:
            tunits = 'days since 1950-01-01 00:00:00 UTC'
        
        fname = f'{outpath}/{outname}_EOBS_{method}.nc'
        if self.info['for_swb2']: ds = nc.Dataset(fname, 'w', format = "NETCDF3_CLASSIC")
        else: ds = nc.Dataset(fname, 'w')
        
        ## General metadata
        ds.description = description
        ds.source = "E-OBS v24.0"
        ds.start_day = start_day
        # ds.author = "paolocolombo1996@gmail.com"
        ds.reference_system = "WGS84"
        ds.proj4_string = "+proj=lonlat +datum=WGS84 +no_defs"
        
        ## Dimensions
        ds.createDimension('x', len(res['lo']))
        ds.createDimension('y', len(res['la']))
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
    
    def save_arcgrid(self):
        self.paths['outpath']
    
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
        if (method == 'cut_area') | (method == 'cut_areatime'): 
            return self.netcdf['latitude'][idx_lat]
        else:
            return self.netcdf['latitude'][:]
        
    def get_lon(self, method, idx_lon = None):
        if (method == 'cut_area') | (method == 'cut_areatime'): 
            return self.netcdf['longitude'][idx_lon]
        else:
            return self.netcdf['longitude'][:]
    
    def get_var(self, method, idx_time = None, idx_lat = None, idx_lon = None):
        if method == 'raw':
            df = np.ma.getdata(self.netcdf[self.info['var']][:, :, :])         
        elif method == 'cut_time':
            df = np.ma.getdata(self.netcdf[self.info['var']][idx_time, :, :])
        elif method == 'cut_area':
            df = np.ma.getdata(self.netcdf[self.info['var']][:, idx_lat, idx_lon])
        elif method == 'cut_areatime':
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
        
    # def transf_daymettime(self):
        
        