# -*- coding: utf-8 -*-
"""
EOBSobject class definition

@author: paolo
"""

import glob
import netCDF4 as nc

class EOBSobject():
    
    def __init__(self, inpath, var, outname = 'none', folder = True):
        #Store the info
        self.info = { 'var': var }
        self.setoutname(outname)
        #Store the path
        if folder: self.paths = { 'inpath': self.find_path(inpath, var) }
        else: self.paths = { 'inpath': inpath }
        #Store the units
        self.units = {
            # 'variable': #save here the unit of the variable
            }
    
    def load(self):
        self.netcdf = nc.Dataset(self.path)
        
    #---------------------------------------------------------
    #NETCDF section
    
    def cut_area(self, coord, save = True, internal = False):
        #self.netcdf
        #df = self.netcdf[self.info['var']][:, :, :]
        if save: self.save()
        if internal: return 1 #la, lo, idx_lat, idx_lon
    
    def cut_time(self, start, end, save = True, internal = False):
        self.netcdf
        #get the time of E-OBS from its metadata
        
        if save: self.save()
        #return the variable?
    
    def cut_areatime(self, save = False):
        self.cut_area(save = False, internal = True)
        self.cut_time(save = False, internal = True)
        
        if save: self.save()
    
    def save(self):
        self.path['outpath']
        #here define the variables and dimensions of the output netcdf
    
    #----------------------------------------------------------
    #ASCII section
    
    #anzi, basta definire bene le funzioni cut in modo che ritornino qualcosa
    #utilizzabile sia da una funzione save_netcdf sia da una save_arcgrid
    
    #----------------------------------------------------------
    #General operations
    
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
        
    def set_timeunit(self, units = 'none'):
        if units == 'none':
            units = self.netcdf['time']['units'] #check if this works
        self.units['time'] = units
        
        
        