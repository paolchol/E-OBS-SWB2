# -*- coding: utf-8 -*-
"""
RechargeCalc class definition
Used to calculate the total recharge
It can be called using:
        from Python.RechargeCalc import RechargeCalc

The working directory has to be set in ./E-OBS-SWB2 for this to work

@author: paolo
"""

# import os
import numpy as np
import pandas as pd
import glob
import time
import geopandas as gp

#The directory has to be set in ./E-OBS-SWB2 for this to work
from Python.SWB2output import SWB2output
from Python.custom_functions import repeat_list

class RechargeCalc():
    
    def __init__(self, swb2path, inputpath, sy, ey, cell_area, uniqueid):
        #Initialize the class
        #Provide the paths to:
        # - the swb2 'net_infiltration' netCDF4 file (swb2path)
        # - the folder where the needed inputh files are saved
        #sy: initial year
        #ey: final year
        
        self.paths = {
            "swb2_output": swb2path,
            "input_folder": inputpath
            }
        self.info = { 
            "start_year": sy,
            "end_year": ey,
            "cell_area_m2": cell_area,
            "id": uniqueid
            }
        self.recharges = {}
        self.conditions = {}
    
    def loadinputfiles(self, irr = True, urb = True):
        #Load the files needed
        print('Loading the input files')
        inpath = self.paths['input_folder']
        fls = glob.glob(f'{inpath}/*.csv')
        names = ['indicatori', 'ricarica_irrigua', 'ricarica_urbana']
        k = []
        for name in names:
            for i in range(len(fls)):
                if(fls[i].find(name) != -1):
                    print(f'{name} file found')
                    k += [i]
        #Get the main indicator file and insert the indicator column
        ind = pd.read_csv(fls[k[0]])
        ind = self.insertind(ind, ind['row'], ind['column'])
        #Store the input files inside the object
        self.input = {
            'ind': ind
            }
        if (irr): self.input['irr'] = pd.read_csv(fls[k[1]])
        if (urb): self.input['urb'] = pd.read_csv(fls[k[2]])
        self.conditions['irr'] = irr
        self.conditions['urb'] = urb
    
    def meteoricR(self, SPs, export = False):
        #Compute the meteoric recharge dataframe
        #Provide the stress periods definition (SPs)
        print('Meteoric recharge dataframe creation')
        start = time.time()
        f = SWB2output(self.paths['swb2_output'])
        rmeteo3d = f.SP_sum(SPs, units = 'ms') #return the SP sum directly in m/s
        f.close()
        
        rmeteo = pd.DataFrame(rmeteo3d[0, :, :])
        rmeteo.insert(0, 'nrow', rmeteo.index.values)
        rmeteo = pd.melt(rmeteo, id_vars = 'nrow', var_name = 'ncol',
         value_name = 'SP1')
        rmeteo['nrow'] = rmeteo['nrow'] + 1
        rmeteo['ncol'] = rmeteo['ncol'] + 1
        rmeteo = self.insertind(rmeteo, rmeteo['nrow'], rmeteo['ncol'])
        
        for i in range(1, rmeteo3d.shape[0]):
            df = pd.DataFrame(rmeteo3d[i, :, :])
            df.insert(0, 'nrow', df.index.values)
            df = pd.melt(df, id_vars = 'nrow', var_name = 'ncol',
             value_name = f'SP{i+1}')
            if f'SP{i+1}' not in rmeteo.columns:
                rmeteo.insert(len(rmeteo.columns), f'SP{i+1}', df[f'SP{i+1}'])
        
        #Save the variables
        self.info['SPs'] = SPs
        self.info['nSP'] = rmeteo3d.shape[0] #number of stress periods
        self.recharges['rmeteo'] = rmeteo
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export:
            outpath = self.paths['outpath'] if 'outpath' in self.paths else self.paths['input_folder']
            rmeteo.to_csv(f'{outpath}/rmeteo.csv')

    def irrigationR(self, Is, coeffs, specialpath = 'none', export = False):
        #Compute the irrigation recharge dataframe
        #Input data: l/s
        print('Irrigation recharge dataframe creation')
        start = time.time()
        nrep = self.info['end_year'] - self.info['start_year'] + 1
        Is = repeat_list(Is, nrep, True)
        irr = self.input['irr']
        
        #Calculate the irrigated area
        if 'area' not in irr.columns:
            irr.insert(len(irr.columns),'area', 0)
        for i, distr in enumerate(irr['distretto'], 0):
            cond = (self.input['ind']['distretto'] == distr) & (self.input['ind']['zona_agricola'] == 1)
            area = sum(cond) * self.info['cell_area_m2']
            irr.loc[i, 'area'] = area
        
        #Calculate the discharge Q in m/s
        if 'Q_ms' not in irr.columns:
            irr.insert(len(irr.columns),'Q_ms', 0)
        irr['Q_ms'] = irr['portata_concessa_ls'] * 0.001 / irr['area']
        
        #Calculate the irrigation recharge and assign it to each cell
        rirr = self.input['ind'].loc[:,(self.info['id'], 'distretto', 'zona_agricola')]
        for i, I in enumerate(Is, 0):
          for j, distr in enumerate(irr['distretto'], 0):
            sp = irr.loc[j, 'speciale'] #'special' code
            if (sp != 1):
                Q = irr['Q_ms'][j]
                cond = (rirr['distretto'] == distr) & (rirr['zona_agricola'] == 1)
                if f'SP{i+1}' not in rirr.columns:
                    rirr.insert(len(rirr.columns), f'SP{i+1}', 0)
                if (sp != 2):
                    rirr.loc[cond, f'SP{i+1}'] = Q * I * coeffs['1'] * coeffs['2']
                else:
                    rirr.loc[cond, f'SP{i+1}'] = Q * I * coeffs['1'] * coeffs['2'] * coeffs['3']
        
        #Assign the provided "special" recharge to the "special" districts
        if (specialpath != 'none'):
            sp_rirr = pd.read_csv(specialpath)
            sdistr = irr.loc[irr['speciale'] == 1, 'distretto']
            for s in sdistr:
                cond = (rirr['distretto'] == s) & (rirr['zona_agricola'] == 1)
                rirrcol = rirr.columns[3:]
                spcol = sp_rirr.columns[1:]
                rirr.loc[cond, rirrcol] = sp_rirr.loc[sp_rirr['distretto'] == s, spcol].values
        
        #Save the variables
        self.recharges['rirr'] = rirr
        self.paths['special_irr'] = specialpath
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export:
            outpath = self.paths['outpath'] if 'outpath' in self.paths else self.paths['input_folder']
            rirr.to_csv(f'{outpath}/rirr.csv')

    def urbanR(self, coeff_urb, export = False):
        #Compute the urban recharge dataframe
        print('Urban recharge dataframe creation')
        start = time.time()
        
        urb = self.input['urb']
        #Calculate the urban area
        if 'area' not in urb.columns:
            urb.insert(1,'area', 0)
        for i, com in enumerate(urb['nome_com'], 0):
            cond = (self.input['ind']['nome_com'] == com) & (self.input['ind']['zona_urbana'] == 1)
            area = sum(cond) * self.info['cell_area_m2']
            urb.loc[i, 'area'] = area
           
        #Calculate the urban recharge and assign it to each cell
        rurb = self.input['ind'].loc[:, (self.info['id'], 'nome_com', 'zona_urbana')]
        for i, colname in enumerate(urb.columns[2:], 0):
            for com in urb['nome_com']:
                if f'SP{i+1}' not in rurb.columns:
                    rurb.insert(len(rurb.columns), f'SP{i+1}', 0)
                Q = urb.loc[urb['nome_com'] == com, colname].values.item()
                A = urb.loc[urb['nome_com'] == com, 'area'].values.item()
                cond = (rurb['nome_com'] == com) & (rurb['zona_urbana'] == 1)
                rurb.loc[cond, f'SP{i+1}'] = (Q / A) * coeff_urb
        
        #Save the variables
        self.recharges['rurb'] = rurb
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export:
            outpath = self.paths['outpath'] if 'outpath' in self.paths else self.paths['input_folder']
            rurb.to_csv(f'{outpath}/rurb.csv')

    def totalR(self, meteopar = None, irrpar = None, urbpar = None, export = False):
        #Sum the recharge components
        print('Total recharge dataframe creation')
        start = time.time()
        #Check if the partial recharges are already computed
        keys = ['rmeteo']
        if 'rmeteo' not in self.recharges:
            self.meteoricR(meteopar['SPs'])
        if self.conditions['irr']:
            keys += ['rirr']
            if 'rirr' not in self.recharges:
                self.irrigationR(irrpar['Is'], irrpar['coeffs'], irrpar['spath'])
        if self.conditions['urb']:
            keys += ['rurb']
            if 'rurb' not in self.recharges:
                self.urbanR(urbpar['coeff_urb'])
        #Sum the recharges
        tool = self.input['ind'].loc[:, self.info['id']]
        tool3d = np.zeros((len(keys), len(tool), self.info['nSP']))
        for i, k in enumerate(keys):
            loc = self.findSPcol(self.recharges[k].columns, self.info['id'])
            toolr = pd.merge(tool, self.recharges[k].loc[:, loc],
                             how = 'left', on = self.info['id'])
            tool3d[i, :, :] = toolr.iloc[:, 1:]
        toolsum = pd.DataFrame(np.sum(tool3d, axis = 0))
        toolsum.columns = toolr.columns[1:]
        toolsum[self.info['id']] = tool
        rtot = self.input['ind'].loc[:, ('row', 'column', self.info['id'])]
        rtot = pd.merge(rtot, toolsum, how = 'left', on = self.info['id'])
        
        #Store the variable
        self.recharges['rtot'] = rtot
        end = time.time()
        print(f'Elapsed time: {round(end - start, 2)} s')
        if export:
            outpath = self.paths['outpath'] if 'outpath' in self.paths else self.paths['input_folder']
            rtot.to_csv(f'{outpath}/rtot.csv')
    
    def findSPcol(self, col, ind):
        names = [ind]
        for name in col:
            if name.find('SP') != -1:
                names += [name]
        return names
    
    def setoutpath(self, outpath):
        if outpath == 'none':
            outpath = self.paths['outpath'] if ('outpath' in self.paths) else self.paths['input_folder']
        else:
            self.paths['outpath'] = outpath
        return outpath
    
    def insertind(self, df, r, c, pos = 0, name = 'none'):
        name = self.info['id'] if name == 'none' else name
        newc = []
        for i in range(len(r)):
            newc += [f'{r[i]}X{c[i]}']
        if (name not in df.columns):
            df.insert(pos, name, newc)
        else:
            df[name] = newc
        return df
    
    def getdf(self, var, tag):
        #Get the df
        if (var == 'input'):
            df = self.input[tag]
        elif (var == 'recharge'):
            df = self.recharges[tag]
        return df

    def export(self, var, tag, outpath = 'none', fileext = 'csv'):
        #Exports the recharges or other data of the class        
        #var: 'input', 'recharge'
        #tag:
        # if var: input
        # - 'ind', 'irr', 'urb'
        # if var: recharge
        # - 'rmeteo', 'rirr', 'rurb', 'rtot'
        #Set the output path
        outpath = self.setoutpath(outpath)
        df = self.getdf(var, tag)
        #Export the file
        df.to_csv(f'{outpath}/{tag}.{fileext}')
    
    def georef(self, var, tag, coordpath, proj = 'none', outpath = 'none'):
        #Export a shapefile of the selected dataframe
        #var: 'input', 'recharge'
        #tag:
        # if var: input
        # - 'ind', 'irr', 'urb'
        # if var: recharge
        # - 'rmeteo', 'rirr', 'rurb', 'rtot'
        #coordpath:
        # path to a .csv file with columns 'row', 'column', self.info['id'], 'x', 'y'
        outpath = self.setoutpath(outpath)
        coord = pd.read_csv(coordpath)
        coord = self.insertind(coord, coord['row'], coord['column'], name = self.info['id'])
        coord = coord.loc[:, (self.info['id'], 'X', 'Y')]
        tool = pd.merge(self.getdf(var, tag), coord, on = self.info['id'])
        geodf = gp.GeoDataFrame(tool, geometry = gp.points_from_xy(tool['X'], tool['Y']))
        geodf.crs = proj if proj != 'none' else '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
        geodf.to_file(f'{outpath}/{tag}.shp', driver = 'ESRI Shapefile')
        print(f'Shapefile saved in {outpath} as {tag}.shp')
        
# def save(self, outpath):
#     #Save a pickle of the class
#prova