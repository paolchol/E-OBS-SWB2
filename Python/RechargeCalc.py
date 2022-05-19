# -*- coding: utf-8 -*-
"""
RechargeCalc class definition
Used to calculate the total recharge
It can be called using:
        from Python.RechargeCalc import RechargeCalc

The working directory has to be set in ./E-OBS-SWB2 for this to work

@author: paolo
"""

import glob
import numpy as np
import pandas as pd
# import sys
import time

#The directory has to be set in ./E-OBS-SWB2 for this to work
from Python.SWB2output import SWB2output

class RechargeCalc():
    
    def __init__(self, swb2path, inputpath, sy, ey, cell_area, uniqueid, nSP):
        """
        Initialize the class
        Provide the paths to:
         - the swb2 'net_infiltration' netCDF4 file (swb2path)
         - the folder where the needed input files are saved (inputpath)
        Other variables needed:
         - sy: initial year
         - ey: final year
         - cell_area: area of the cell in m2
         - uniqueid: name of the unique id column in the "indicatori" file
        """
        
        self.paths = {
            "swb2_output": swb2path,
            "input_folder": inputpath
            }
        self.info = { 
            "start_year": sy,
            "end_year": ey,
            "cell_area_m2": cell_area,
            "id": uniqueid,
            "nSP": nSP
            }
        self.recharges = {}
        self.conditions = {}
    
    def load_inputfiles(self, meteo = True, irr = True, urb = True):
        #Load the files needed
        print('Loading the input files')
        print('-----------------------')
        inpath = self.paths['input_folder']
        fls = glob.glob(f'{inpath}/*.csv')
        names = ['indicatori', 'ricarica_irrigua', 'extractions']
        k = []
        for name in names:
            for i in range(len(fls)):
                if(fls[i].find(name) != -1):
                    print(f'{name} file found')
                    self.paths[name] = fls[i]
                    k += [i]
        #Get the main indicator file and insert the custom indicator column
        ind = pd.read_csv(fls[k[0]])
        ind = self.insert_ind(ind, ind['row'], ind['column'])
        #Store the input files inside the object
        self.input = { 'ind': ind }
        if (irr): self.input['irr'] = pd.read_csv(fls[k[1]])
        if (urb): self.input['urb'] = pd.read_csv(fls[k[2]])
        self.conditions['meteo'] = meteo
        self.conditions['irr'] = irr
        self.conditions['urb'] = urb
    
    #----------------------------------------------------------------------
    #Recharges calculation
    
    def meteoricR(self, SPs, units = 'ms', fixrow = 1, fixcol = 1, export = False, ret = False):
        """
        Compute the meteoric recharge dataframe
        SPs: Stress periods definition, in days
            e.g. SPs = [90, 76, 92, 107] represents 4 SP of length 90, 76..
        """
        print('Meteoric recharge dataframe creation')
        print('------------------------------------')
        start = time.time()
        SPs = self.set_SPs(SPs, 1)
        f = SWB2output(self.paths['swb2_output'])
        rmeteo3d = f.SP_sum(SPs, units = units, retval = True) #return the SP sum directly in m/s
        f.close()
        
        rmeteo = pd.DataFrame(rmeteo3d[0, :, :])
        rmeteo.insert(0, 'nrow', rmeteo.index.values)
        rmeteo = pd.melt(rmeteo, id_vars = 'nrow', var_name = 'ncol', value_name = 'SP1')
        rmeteo = self.insert_ind(rmeteo, rmeteo['nrow'], rmeteo['ncol'], fixrow, fixcol)
        
        for i in range(1, rmeteo3d.shape[0]):
            df = pd.DataFrame(rmeteo3d[i, :, :])
            df.insert(0, 'nrow', df.index.values)
            df = pd.melt(df, id_vars = 'nrow', var_name = 'ncol', value_name = f'SP{i+1}')
            df = self.insert_ind(df, df['nrow'], df['ncol'], fixrow, fixcol)
            if f'SP{i+1}' not in rmeteo.columns:
                rmeteo = rmeteo.join(df.loc[:,[self.info['id'], f'SP{i+1}']].set_index(self.info['id']), on = self.info['id'])
                
        lastSP = i+1
        if lastSP != self.info['nSP']:
            print(f"Replicating columns to create a dataframe of {self.info['nSP']} SPs")           
            x = lastSP
            y = self.info['nSP']
            k, d = divmod(y/x, 1)
            k = round(k)
            tool = rmeteo.set_index('indicatore').copy()
            concatenated = pd.concat([tool[self.find_SPcol(tool)]]*k, axis = 1)
            sps = self.find_SPcol(rmeteo, 'indicatore', True)
            df = rmeteo.loc[:, sps[0:round(x*d)+1]].copy()
            joined = concatenated.join(df.set_index(self.info['id']), on = self.info['id'], rsuffix = '_new')
            cc = [f'SP{i+1}' for i in range(len(joined.columns))]
            joined.columns = cc
            joined = joined.reset_index(level = 0)
            rmeteo = joined.copy()
        
        #Save the variables
        self.info['SPs'] = SPs
        self.recharges['rmeteo'] = rmeteo
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export: self.export('recharge', 'rmeteo')
        if ret: return rmeteo
    
    def irrigationR(self, coeffs, specialpath = 'none', export = False,
                    multicoeff = False, splist = None):
        """
        coeffs: dictionary if multicoeff is False
                pandas.DataFrame if multicoeff is True
        splist: list of SPs in which to use the second line of coefficients
        """
        print('Irrigation recharge dataframe creation')
        print('--------------------------------------')
        start = time.time()
        irr = self.input['irr']
        if specialpath != 'none':
            sp_irr = pd.read_csv(specialpath)
            self.input['special_irr'] = sp_irr
        
        #Calculate the irrigated area
        if 'area' not in irr.columns:
            irr.insert(len(irr.columns),'area', 0)
        for distr in irr['distretto']:
            cond = (self.input['ind']['distretto'] == distr) & (self.input['ind']['zona_agricola'] == 1)
            area = sum(cond) * self.info['cell_area_m2']
            irr.loc[irr['distretto'] == distr, 'area'] = area
        
        #Calculate the irrigation recharge and assign it to each cell
        rirr = self.input['ind'].loc[:, (self.info['id'], 'distretto', 'zona_agricola')]
        # tool = pd.DataFrame(np.zeros(len(rirr.index), len(self.find_SPcol(irr.columns))))
        for sp in self.find_SPcol(irr.columns):
            if sp not in rirr.columns:
                rirr.insert(len(rirr.columns), sp, 0)
                
                # rirr.join(pd.DataFrame({sp: [0 for x in range(len(rirr.index)]}))
                # pd.concat((rirr, pd.DataFrame({sp: [0 for x in range(0, len(rirr.index))]})), axis = 1)
            if multicoeff:
                x = 1 if sp in splist else 0
                RISP = coeffs['RISP'][x]
                P = coeffs['P'][x]
                K = 1 - coeffs['E'][x] - coeffs['R'][x]
            else:
                RISP = coeffs['RISP']
                P = coeffs['P']
                K = 1 - coeffs['E'] - coeffs['R']
            for distr in irr['distretto']:
                code = irr.loc[irr['distretto'] == distr, 'code'].values[0]
                cond = (rirr['distretto'] == distr) & (rirr['zona_agricola'] == 1)
                if code != 1:
                    Q = float(irr.loc[irr['distretto'] == distr, sp].values[0])
                    A = irr.loc[irr['distretto'] == distr, 'area'].values[0]
                    rirr.loc[cond, sp] = (Q * RISP)/(A * P) * K
                else:
                    Q = float(sp_irr.loc[sp_irr['distretto'] == distr, sp].values[0])
                    rirr.loc[cond, sp] = Q
        
        #Store the variables
        self.recharges['rirr'] = rirr.copy() #save a de-fragmented dataframe
        self.paths['special_irr'] = specialpath
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export: self.export('recharge', 'rirr')
    
    def urbanR(self, coeff, col = None, valcol = None, option = None,
               areas = False, export = False):
        """
        Compute the urban recharge dataframe
        The urban recharge is calculated as a fraction of the pumped volumes. It
        is due to the losses from the extraction pumps and pipes.
        
        coeff: coefficient to apply to the extractions
        col: list object, containing strings
        valcol: list object, containing any format compatible with the values
                to check
        option: list object, containing 1 if the condition wanted is AND, or
                0 if the condition wanted is OR
        """
        print('Urban recharge dataframe creation')
        print('---------------------------------')
        start = time.time()
        single_cond = False if len(col) > 1 else True
        if single_cond:
            cond = self.input['ind'][col[0]] == valcol[0]
        else:
            cond = self.input['ind'][col[0]] == valcol[0]
            for i in range(1, len(col)):                
                if option[i-1] == 1:
                    cond = (cond) & (self.input['ind'][col[i]] == valcol[i])
                else:                    
                    cond = (cond) | (self.input['ind'][col[i]] == valcol[i])
        
        urb = self.input['urb'].copy()
        rurb = self.input['ind'].loc[:, [self.info['id'], 'nome_com'] + list(set(col))]
        tool = pd.DataFrame(np.zeros((len(rurb), self.info['nSP'])), index = rurb[self.info['id']])
        tool.columns = [f'SP{i+1}' for i in range(self.info['nSP'])]
        rurb = rurb.join(tool, on = self.info['id'])
        
        for com in urb['nome_com']:
            idx = (cond) & (self.input['ind']['nome_com'] == com)
            A = sum(idx) * self.info['cell_area_m2']
            for sp in self.find_SPcol(urb):
                E = abs(urb.loc[urb['nome_com'] == com, sp].values.item())
                rurb.loc[idx, sp] = E / A * coeff
        
        if areas:
            area = [sum((cond) & (self.input['ind']['nome_com'] == com)) * self.info['cell_area_m2'] for com in urb['nome_com']]
            urb.insert(1, 'area', area)
            self.input['urb'] = urb
        
        #Store the variables
        self.recharges['rurb'] = rurb.copy()
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export: self.export('recharge', 'rurb')
    
    def totalR(self, meteopar = None, irrpar = None, urbpar = None, export = False):
        """
        Sum the recharge components
        """
        print('Total recharge dataframe creation')
        print('---------------------------------')
        start = time.time()
        #Check if the partial recharges are already computed
        keys = []
        if self.conditions['meteo']:
            keys += ['rmeteo']
            if 'rmeteo' not in self.recharges:
                self.meteoricR(meteopar['SPs'])
        if self.conditions['irr']:
            keys += ['rirr']
            if 'rirr' not in self.recharges:
                self.irrigationR(irrpar['coeffs'], irrpar['spath'])
        if self.conditions['urb']:
            keys += ['rurb']
            if 'rurb' not in self.recharges:
                self.urbanR(urbpar['coeff'], urbpar['col'])
        #Sum the recharges
        tool = self.input['ind'].loc[:, self.info['id']]
        tool3d = np.zeros((len(keys), len(tool), self.info['nSP']))
        for i, k in enumerate(keys):
            loc = self.find_SPcol(self.recharges[k].columns, self.info['id'], True)
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
        if export: self.export('recharge', 'rtot')
    
    #-----------------------------------------------------------------------
    #Operations on the recharges
    
    # def modify_recharge(self, var, tag, cond, operation):
    #     """
    #     Mock-up of a possible function to operate directly on the recharges
    #     """
    #     df = self.get_df(var, tag)
    #     df[cond] = operation(df[cond])
    #     self.recharges[tag] = df
    
    def modify_recharge(self, var, tag, coeff, single_cond = True,
                        multi_cond = False, col = None, valcol = None):
        """
        Modifies the values of the cells that have col == valcol
        multiplying them by a coefficient (coeff)
        
        col and valcol can be provided as lists, by setting single_cond to False
        and multi_cond as True
        """
        if single_cond:
            cond = self.input['ind'][col] == valcol
        elif multi_cond:
            cond = self.input['ind'][col[0]] == valcol[0]
            for i in range(1, len(col)):
                cond = (cond) & (self.input['ind'][col[i]] == valcol[i])
        df = self.get_df(var, tag)
        idx = self.input['ind'].loc[cond, self.info['id']]
        idx2 = df[self.info['id']].isin(idx)
        df.loc[idx2, self.find_SPcol(df)] = df.loc[idx2, self.find_SPcol(df)] * coeff
        self.recharges[tag] = df
        
    def add_attibute():
        #select a recharge
        #add attributes from input['ind'] to the selected recharge df
        pass
    
    #-----------------------------------------------------------------------
    #General functions
    
    def double_cond(self, c1, c2):
        cond = c1
        for i in range(len(c2)):
            cond = c1 & c2[i]
        return cond
    
    def export(self, var, tag, fileext = 'csv', outpath = 'none', outname = 'none',
               withcoord = False, coordpath = None):
        """
        Exports the recharges or other data of the class        
        var: 'input', 'recharge'
        tag:
         if var = input
         - 'ind', 'irr', 'urb'
         if var = recharge
         - 'rmeteo', 'rirr', 'rurb', 'rtot'
        fileext: file extention wanted. Default: csv
        outpath: path to a wanted output folder. Default: variable 'outpath'
         defined previously
        """
        start = time.time()
        #Set the output path
        outpath = self.set_outpath(outpath)
        df = self.get_df(var, tag)
        if withcoord: 
            coord = pd.read_csv(coordpath)
            coord = self.insert_ind(coord, coord['row'], coord['column'], name = self.info['id'])
            coord = coord.loc[:, (self.info['id'], 'X', 'Y')]
            df = pd.merge(coord, df, on = self.info['id'])
        #Export the file
        outname = tag if outname == 'none' else outname
        df.to_csv(f'{outpath}/{outname}.{fileext}', index = False)
        end = time.time()
        print(f'{end-start} s')
        
    def find_SPcol(self, col, ind = None, indname = False):
        names = [ind] if indname else []
        for name in col:
            if name.find('SP') != -1:
                names += [name]
        return names
    
    def get_df(self, var, tag):
        #Get the df
        if var == 'input':
            df = self.input[tag]
        elif var == 'recharge':
            df = self.recharges[tag]
        return df
    
    def georef(self, var, tag, coordpath, crs = 'epsg:4326', proj = 'none', outpath = 'none',
               outname = 'none', dropcoord = False):
        """
        Export a shapefile of the selected dataframe
        var (str): 'input', 'recharge'
        tag (str):
         if var: input
         - 'ind', 'irr', 'urb'
         if var: recharge
         - 'rmeteo', 'rirr', 'rurb', 'rtot'
        coordpath:
         path to a .csv file with columns 'row', 'column', self.info['id'], 'x', 'y'
        outpath: path to a wanted output folder. Default: variable 'outpath'
         defined previously
        """
        #if 'geopandas' not in sys.modules: 
        import geopandas as gp
        
        start = time.time()
        outpath = self.set_outpath(outpath)
        outname = outname if outname != 'none' else f'{tag}'
        coord = pd.read_csv(coordpath)
        coord = self.insert_ind(coord, coord['row'], coord['column'], name = self.info['id'])
        coord = coord.loc[:, (self.info['id'], 'X', 'Y')]
        tool = pd.merge(coord, self.get_df(var, tag), on = self.info['id'])
        # newtool = tool.copy()
        geodf = gp.GeoDataFrame(tool, geometry = gp.points_from_xy(tool['X'], tool['Y']))
        # geodf.crs = proj if proj != 'none' else '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
        geodf.set_crs(crs, inplace = True)
        geodf.index = geodf[self.info['id']]
        geodf.drop(self.info['id'], axis = 1, inplace = True)
        if dropcoord: geodf.drop(['X', 'Y'], axis = 1, inplace = True)

        #This operation (to_file) takes too much time, try to address the issues that causes it and
        #fix them
        geodf.to_file(f'{outpath}/{outname}.shp')
        end = time.time()
        print(f'Shapefile saved in {outpath} as {outname}.shp')
        print(f'Elapsed time: {round(end-start, 2)} s')
    
    def insert_ind(self, df, r, c, fixrow = 0, fixcol = 0, pos = 0, name = 'none'):
        name = self.info['id'] if name == 'none' else name
        if (fixrow != 0) | (fixcol != 0):
            r = r + fixrow
            c = c + fixcol
        newc = []
        for i in range(len(r)):
            newc += [f'{r[i]}X{c[i]}']
        if (name not in df.columns):
            df.insert(pos, name, newc)
        else:
            df[name] = newc
        return df
    
    def set_SPs(self, SPs, c = None):
        """
        Sets the stress periods duration as a variable of the class
        If c = 1 returns the cumulative sum of the stress periods duration
        SPs (list):
            A list containing the length in days of the stress periods inside
            one year
        """
        self.info['SPs'] = SPs
        if c == 1: return np.cumsum(SPs)
    
    def set_tags(self, district = None, agr_zone = None, urb_zone = None):
        pass
    
    def set_outpath(self, outpath):
        if outpath == 'none':
            outpath = self.paths['outpath'] if 'outpath' in self.paths else self.paths['input_folder']
        else:
            self.paths['outpath'] = outpath
        return outpath
    
    def save(self, outpath):
         #Save a pickle of the object
         pass