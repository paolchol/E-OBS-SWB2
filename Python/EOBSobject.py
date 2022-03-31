# -*- coding: utf-8 -*-
"""
EOBSobject class definition

@author: paolo
"""

import glob
import netCDF4 as nc
import numpy as np

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
        df = self.netcdf[self.info['var']][:, idx_lat, idx_lon]
        df = np.ma.getdata(df)
        res = {
            'df': df,
            'la': la[idx_lat],
            'lo': lo[idx_lon]
            }
        if internal: return res
        if save: self.save_netcdf(res, 'cut_area')
    
    def cut_time(self, start, end, save = True, internal = False):
        # self.netcdf
        #get the time of E-OBS from its metadata
        
        #inserire condizione sul metodo di registrazione del tempo
        #time = eobs_todaymet(np.where(tyR == year)[0])
        
        #opzione: file unico ritagliato negli anni voluti
        #file singoli anno per anno
        
        #serve qui un ciclo for sul tempo
        #che a ogni passaggio chiami save_netcdf
        
        # res = {
        #     'df': df,
        #     'la': la,
        #     'lo': lo,
        #     'year': year,
        #     'time': time
        #     }
        
        if save: self.save_netcdf()
        #return the variable?
    
    def cut_areatime(self, save = False):
        self.cut_area(save = False, internal = True)
        self.cut_time(save = False, internal = True)
        
        if save: self.save_netcdf()
    
    def save_netcdf(self, res, method = 'none'):
        outpath = self.path['outpath']
        outname = self.info['outname']
        
        if self.info['for_swb2']:
            res['la'][::-1].sort()
            res['df'] = np.flip(res['df'], axis = 1)
            res['df'] = np.around(res['df'], 1)
        
        if method == 'none':
            tag = 'raw'
            #save the same
            pass
        
        if method == 'cut_area':
            tag = method
            description = "Clip in space of the E-OBS dataset"
            # time = estrai tempo dal netcdf
            
        if method == 'cut_time':
            description = "Clip in time of the E-OBS dataset"
            tag = method
            start_day = f"01/01/{res['year']}"
            # end_day = calcola l'end day
            time = res['time']
        
        if method == 'cut_areatime':
            tag = method
            
        fname = f'{outpath}/{outname}_EOBS_{tag}.nc'
        ds = nc.Dataset(fname, 'w', format = "NETCDF3_CLASSIC")
        
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
        time.units = 'days since 1980-01-01 00:00:00 UTC' #try to keep the same numeration as E-OBS
        #Day of the year
        yearday = ds.createVariable('yearday', 'h', ('time'))
        #Value
        value = ds.createVariable(outname, 'f4', ('time','y','x'), fill_value = -9999)
        value.units = self.info['units']
        value.missing_value = -9999.0
        value.coordinates = 'lat lon'
        
        ## Fill the variables
        x[:] = res['lo']
        y[:] = res['la']
        tout = time
        time[:] = tout
        yearday[:] = range(1, 366)
        value[:] = res['df']
        
        #Close the file
        ds.close()
    
    def save_arcgrid(self):
        self.path['outpath']
    
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
    
    def set_outname(self, outname):
        self.info['outname'] = outname if outname != 'none' else self.info['var']
    
    def set_outpath(self, outpath, folder):
        tool = self.paths['inpath']
        if not folder:
            tool = tool.split('/')[:-1]
            tool = '/'.join(tool)
        else: tool = tool.split('\\')[0]
        self.path['outpath'] = outpath if outpath != 'none' else tool
    
    def set_timeunit(self, units = 'none'):
        if units == 'none':
            units = self.netcdf['time']['units'] #check if this works
        self.units['time'] = units
        
        
        