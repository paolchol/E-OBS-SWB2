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
    
    def load_inputfiles(self, irr = True, urb = True):
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
                    self.paths[name] = fls[i]
                    k += [i]
        #Get the main indicator file and insert the custom indicator column
        ind = pd.read_csv(fls[k[0]])
        ind = self.insert_ind(ind, ind['row'], ind['column'])
        #Store the input files inside the object
        self.input = {
            'ind': ind
            }
        if (irr): self.input['irr'] = pd.read_csv(fls[k[1]])
        if (urb): self.input['urb'] = pd.read_csv(fls[k[2]])
        self.conditions['irr'] = irr
        self.conditions['urb'] = urb
    
    #----------------------------------------------------------------------
    #Recharges calculation
    
    def meteoricR(self, SPs, export = False):
        #Compute the meteoric recharge dataframe
        #Provide the stress periods definition (SPs)
        print('Meteoric recharge dataframe creation')
        start = time.time()
        SPs = self.set_SPs(SPs, 1)
        f = SWB2output(self.paths['swb2_output'])
        rmeteo3d = f.SP_sum(SPs, units = 'ms', retval = True) #return the SP sum directly in m/s
        f.close()
        
        return rmeteo3d
        
        rmeteo = pd.DataFrame(rmeteo3d[0, :, :])
        rmeteo.insert(0, 'nrow', rmeteo.index.values)
        rmeteo = pd.melt(rmeteo, id_vars = 'nrow', var_name = 'ncol', value_name = 'SP1')
        rmeteo['nrow'] = rmeteo['nrow'] + 1
        rmeteo['ncol'] = rmeteo['ncol'] + 1
        rmeteo = self.insert_ind(rmeteo, rmeteo['nrow'], rmeteo['ncol'])
        
        for i in range(1, rmeteo3d.shape[0]):
            df = pd.DataFrame(rmeteo3d[i, :, :])
            df.insert(0, 'nrow', df.index.values)
            df = pd.melt(df, id_vars = 'nrow', var_name = 'ncol', value_name = f'SP{i+1}')
            # df = self.insert_ind(df, 'nrow', 'ncol')
            if f'SP{i+1}' not in rmeteo.columns:
                rmeteo.insert(len(rmeteo.columns), f'SP{i+1}', df[f'SP{i+1}'])
                #rmeteo.join(pd.DataFrame({f'SP{i+1}': df[f'SP{i+1}']}), on = self.info['ind'])
        
        lastSP = i+1
        if lastSP != self.info['nSP']:
            sps = self.find_SPcol(rmeteo.columns)
            k = 0
            for i in range(lastSP, self.info['nSP']):
                # rmeteo.join(pd.DataFrame({f'SP{i+1}': rmeteo.loc[:, sps[k]]}))
                rmeteo.insert(len(rmeteo.columns), f'SP{i+1}', rmeteo.loc[:, sps[k]])
                k = k+1 if k < len(sps)-1 else 0
        
        #Save the variables
        self.info['SPs'] = SPs
        # self.info['nSP'] = rmeteo3d.shape[0] #number of stress periods
        self.recharges['rmeteo'] = rmeteo
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export: self.export('recharge', 'rmeteo')
    
    def irrigationR(self, coeffs, specialpath = 'none', export = False,
                    multicoeff = False, splist = None):
        """
        coeffs: dictionary if multicoeff is False, pandas.DataFrame if multicoeff
                is True
        splist: list of SPs in which to use the second line of coefficients
        """
        print('Irrigation recharge dataframe creation')
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
        for sp in self.find_SPcol(irr.columns):
            if sp not in rirr.columns:
                rirr.insert(len(rirr.columns), sp, 0)
                # rirr.join(pd.DataFrame({sp: [0 for x in range(0, len(rirr.index))]}))
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
    
    def urbanR(self, coeff_urb, export = False):
        """
        Compute the urban recharge dataframe
        The urban recharge is calculated as a fraction of a provided discharge
        
        """
        print('Urban recharge dataframe creation')
        start = time.time()
        
        urb = self.input['urb']
        #Calculate the urban area
        if 'area' not in urb.columns:
            urb.insert(1,'area', 0)
        for i, com in enumerate(urb['nome_com'], 0):
            cond = self.double_cond(self.input['ind']['nome_com'] == com, self.input['ind']['zona_urbana'] == 1)
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
                cond = self.double_cond(rurb['nome_com'] == com, rurb['zona_urbana'] == 1)
                rurb.loc[cond, f'SP{i+1}'] = (Q / A) * coeff_urb
        
        #Save the variables
        self.recharges['rurb'] = rurb
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export: self.export('recharge', 'rurb')
    
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
                self.irrigationR(irrpar['coeffs'], irrpar['spath'])
        if self.conditions['urb']:
            keys += ['rurb']
            if 'rurb' not in self.recharges:
                self.urbanR(urbpar['coeff_urb'])
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
    
    def insert_ind(self, df, r, c, pos = 0, name = 'none'):
        name = self.info['id'] if name == 'none' else name
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
    
    #-------------------------------------------------------
    #Old functions section
    
    def old_irrigationR(self, Is, coeffs, specialpath = 'none', export = False):
        #Compute the irrigation recharge dataframe
        #Input data: l/s
        print('Irrigation recharge dataframe creation')
        from Python.custom_functions import repeat_list
        start = time.time()
        nrep = self.info['end_year'] - self.info['start_year'] + 1
        Is = repeat_list(Is, nrep, True)
        irr = self.input['irr']
        
        #Calculate the irrigated area
        if 'area' not in irr.columns:
            irr.insert(len(irr.columns),'area', 0)
        for i, distr in enumerate(irr['distretto'], 0):
            cond = self.double_cond(self.input['ind']['distretto'] == distr, self.input['ind']['zona_agricola'] == 1)
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
                cond = self.double_cond(rirr['distretto'] == distr, rirr['zona_agricola'] == 1)
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
                cond = self.double_cond(rirr['distretto'] == s, rirr['zona_agricola'] == 1)
                rirrcol = rirr.columns[3:]
                spcol = sp_rirr.columns[1:]
                rirr.loc[cond, rirrcol] = sp_rirr.loc[sp_rirr['distretto'] == s, spcol].values
        
        #Save the variables
        self.recharges['rirr'] = rirr
        self.paths['special_irr'] = specialpath
        if (specialpath != 'none'): self.input['spirr'] = sp_rirr
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export:
            outpath = self.paths['outpath'] if 'outpath' in self.paths else self.paths['input_folder']
            rirr.to_csv(f'{outpath}/rirr.csv')
        
# def save(self, outpath):
#     #Save a pickle of the object