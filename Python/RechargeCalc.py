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
import os
import pandas as pd
import sys
import time

sys.path.append(os.getcwd())

# Import the SWB2output class
from SWB2output import SWB2output

class RechargeCalc():
    
    def __init__(self, sy, ey, cell_area, uniqueid, nSP, customid = False,
                 meteo = True, irr = True, urb = True):
        """
        Initialize the class. Creates the dictionaries "info", "recharges" and 
        "conditions". Select which recharge components to consider.

        Parameters:
        ----------
        sy : int
            Initial year
        ey : int
            Final year
        cell_area : float
            Area of the cell in square meters
        uniqueid : str
            Name of the unique id column in the "indicatori" file
        nSP : int
            Number of stress periods
        customid : bool
            If you file already has a custom index to address each cell and you
            want to keep it to use it later set this to True.
        meteo : bool
            Consideration of the meteoric recharge component. The default is True.
        irr : bool
            Consideration of the irrigation recharge component. The default is True.
        urb : bool
            DESCRIPTION. The default is True.
        """
        self.info = {
            "start_year": sy,
            "end_year": ey,
            "cell_area_m2": cell_area,
            "id": uniqueid,
            "nSP": nSP,
            "customid": customid
            }
        self.recharges = {}
        self.conditions = {}
        self.sel_recharge(meteo, irr, urb)
    
    def load_inputfiles(self, swb2path, inputpath = None):
        """
        Load the input files needed

        Parameters
        ----------
        swb2path : str
            Path to SWB2's output 'net_infiltration' NetCDF file.
        inputpath : str
            Path to the folder where all the other needed files are stored.
            These files need to be in the same folder. The default is None.
        """
        #Save the paths
        self.paths = {
            "swb2_output": swb2path,
            "input_folder": inputpath
            }
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
        if not self.info['customid']: ind = self.insert_ind(ind, ind['row'], ind['column'])
        #Store the input files inside the object
        self.input = { 'ind': ind }
        if self.conditions['irr']: self.input['irr'] = pd.read_csv(fls[k[1]])
        if self.conditions['urb']: self.input['urb'] = pd.read_csv(fls[k[2]])
    
    #Recharges calculation
    #---------------------
    
    def meteoricR(self, SPs, units = 'ms', fixrow = 1, fixcol = 1, export = False, ret = False):
        """
        Compute the meteoric recharge dataframe
        SPs: list of int
            Stress periods definition, in days
            e.g. SPs = [90, 76, 92, 107] represents 4 SP of length 90, 76.
        fixrow: int, optional
            Row location
        fixcol: int, optional
            Column location
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
        rmeteo['ncol'] = pd.to_numeric(rmeteo['ncol'])
        rmeteo = self.insert_ind(rmeteo, rmeteo['nrow'], rmeteo['ncol'], fixrow, fixcol)
        
        for i in range(1, rmeteo3d.shape[0]):
            df = pd.DataFrame(rmeteo3d[i, :, :])
            df.insert(0, 'nrow', df.index.values)
            df = pd.melt(df, id_vars = 'nrow', var_name = 'ncol', value_name = f'SP{i+1}')
            df = self.insert_ind(df, df['nrow'], df['ncol'], fixrow, fixcol)
            if f'SP{i+1}' not in rmeteo.columns:
                rmeteo = rmeteo.join(df.loc[:,[self.info['id'], f'SP{i+1}']].set_index(self.info['id']), on = self.info['id'])
        lastSP = i+1
        # Replicate columns if needed
        if lastSP < self.info['nSP']:         
            x = lastSP
            rmeteo = self.replicate_columns(x, rmeteo)        
        #Save the variables
        self.info['SPs'] = SPs
        self.recharges['rmeteo'] = rmeteo
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export: self.export('recharge', 'rmeteo')
        if ret: return rmeteo
    
    def irrigationR(self, coeffs, specialpath = 'none', multicoeff = False,
                    splist = None, areas = False, export = False, ret = False):
        """
        coeffs: dictionary
                Contains the coefficients needed to calculate the irrigation
                recharge from the provided discharge in each SP
        splist: list of SPs in which to use the second line of coefficients
        """
        print('Irrigation recharge dataframe creation')
        print('--------------------------------------')
        start = time.time()
        #Load input and parameters
        irr = self.input['irr'].copy()
        if specialpath != 'none':
            sp_irr = pd.read_csv(specialpath)
            self.input['special_irr'] = sp_irr
        if not multicoeff:    
            RISP = coeffs['RISP']
            P = coeffs['P']
            K = 1 - coeffs['E'] - coeffs['R']
        #Create the output dataframe
        cond = self.input['ind']['zona_agricola'] == 1
        rirr = self.input['ind'].loc[:, [self.info['id'], 'distretto', 'zona_agricola']]
        tool = pd.DataFrame(np.zeros((len(rirr), len(self.find_SPcol(irr)))), index = rirr[self.info['id']])
        tool.columns = [f'SP{i+1}' for i in range(len(self.find_SPcol(irr)))]
        rirr = rirr.join(tool, on = self.info['id'])
        #Calculate the irrigation recharge for each district and SP
        for distr in irr['distretto']:
            code = irr.loc[irr['distretto'] == distr, 'code'].values[0]
            idx = (cond) & (self.input['ind']['distretto'] == distr)
            if code != 1: A = sum(idx) * self.info['cell_area_m2']
            for sp in self.find_SPcol(irr):
                if multicoeff:
                    x = 1 if sp in splist else 0
                    RISP = coeffs['RISP'][x]
                    P = coeffs['P'][x]
                    K = 1 - coeffs['E'][x] - coeffs['R'][x]
                if code != 1:
                    Q = float(irr.loc[irr['distretto'] == distr, sp].values[0])
                    rirr.loc[idx, sp] = (Q * RISP)/(A * P) * K
                else:
                    Q = float(sp_irr.loc[sp_irr['distretto'] == distr, sp].values[0])
                    rirr.loc[idx, sp] = Q
        # Replicate columns if needed
        if len(self.find_SPcol(irr))<self.info['nSP']:
            x = len(self.find_SPcol(irr))
            rirr = self.replicate_columns(x, rirr)
        #Store the variables
        if areas:
            area = [sum((cond) & (self.input['ind']['distretto'] == distr)) * self.info['cell_area_m2'] for distr in irr['distretto']]
            irr.insert(1, 'area', area)
            self.input['irr'] = irr
        self.recharges['rirr'] = rirr.copy()
        self.paths['special_irr'] = specialpath
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export: self.export('recharge', 'rirr')
        if ret: return rirr
    
    def urbanR(self, coeff, col = None, valcol = None, option = None,
               areas = False, export = False, ret = False):
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
        tool = pd.DataFrame(np.zeros((len(rurb), len(self.find_SPcol(urb)))), index = rurb[self.info['id']])
        tool.columns = [f'SP{i+1}' for i in range(len(self.find_SPcol(urb)))]
        rurb = rurb.join(tool, on = self.info['id'])
        
        for com in urb['nome_com']:
            idx = (cond) & (self.input['ind']['nome_com'] == com)
            A = sum(idx) * self.info['cell_area_m2']
            for sp in self.find_SPcol(urb):
                E = abs(urb.loc[urb['nome_com'] == com, sp].values.item())
                rurb.loc[idx, sp] = E / A * coeff
        # Replicate columns if needed
        if len(self.find_SPcol(urb))<self.info['nSP']:
            x = len(self.find_SPcol(urb))
            rurb = self.replicate_columns(x, rurb)
        #Store the variables
        if areas:
            area = [sum((cond) & (self.input['ind']['nome_com'] == com)) * self.info['cell_area_m2'] for com in urb['nome_com']]
            urb.insert(1, 'area', area)
            self.input['urb'] = urb
        self.recharges['rurb'] = rurb.copy()
        end = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export: self.export('recharge', 'rurb')
        if ret: return rurb
    
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
            if toolr.shape[1] > tool3d.shape[2]:
                toolr = toolr.loc[:, toolr.columns[:tool3d.shape[2]+1]]
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
    
    #Operations on the recharges
    #---------------------------
    
    def join_external(self, var, tag, extdf, on = 'none', rsuffix = 'none'):
        """
        Function to join external dataframes to the ones stored in the object.
        For example, join the total recharge dataframes of two different periods
        to obtain a single dataframe.
        
        Parameters
        ----------
        var : str
        tag : str
            as defined in get_df
        extdf : pandas.DataFrame object
            the dataframe you want to join to the one extracted from
            RechargeCalc
        """       
        on = self.info['id'] if on == 'none' else on
        df = self.get_df(var, tag)
        if rsuffix != 'none': df = df.join(extdf, on = on, rsuffix = rsuffix)
        else: df = df.join(extdf, on = on)
        return df
    
    def load_component(self, var, tag, df = None, read = False, path = None):
        #Insert the df
        if read: df = pd.read_csv(path)
        if var == 'input':
            self.input[tag] = df
        elif var == 'recharge':
            self.recharges[tag] = df

    def modify_recharge(self, var, tag, coeff, cond = 'pass', single_cond = True,
                        multi_cond = False, col = None, valcol = None):
        """
        Modifies the values of the cells that have col == valcol
        multiplying them by a coefficient (coeff)
        
        col and valcol can be provided as lists, by setting single_cond to False
        and multi_cond as True

        var : str
            as defined in RechargeCalc.get_df
        tag : str
            as defined in RechargeCalc.get_df
        coeff: float
            multiplication coefficient
        cond: str, optional
            if no condition is imposed (i.e. the multiplication coeff
            will be applied to all rows) set this parameter as "null". Default is "pass"
        single_cond: bool, optional
            Set it to True if a single condition on the rows is applied.
            Default is True
        multi_cond: bool, optional
            Set it to True and set single_cond to False if multiple conditions are passed
            Default is False
        col: str of list of str, optional
            the column inside the indicator file where to apply the coeff
        valcol: str or list 
            the values of the respective column inside the indicator
            file where to apply the coeff
        """
        if cond == 'null':
            pass
        elif single_cond:
            cond = self.input['ind'][col] == valcol
        elif multi_cond:
            cond = self.input['ind'][col[0]] == valcol[0]
            for i in range(1, len(col)):
                cond = (cond) & (self.input['ind'][col[i]] == valcol[i])
        df = self.get_df(var, tag)
        if cond != 'null':
            idx = self.input['ind'].loc[cond, self.info['id']]
        else:
            idx = self.input['ind'].loc[:, self.info['id']]
        idx2 = df[self.info['id']].isin(idx)
        df.loc[idx2, self.find_SPcol(df)] = df.loc[idx2, self.find_SPcol(df)] * coeff
        self.recharges[tag] = df
        
    def add_attibute():
        #select a recharge
        #add attributes from input['ind'] to the selected recharge df
        pass
    
    #Export functions
    #----------------
    
    def export(self, var, tag, fileext = 'csv', outpath = None, outname = 'none',
               withcoord = False, coordpath = None):
        """
        Exports the recharges or other data of the class
        
        Parameters
        ----------
        var : str
            As defined in get_df
        tag : str
            As defined in get_df
        fileext : str
            File extention wanted. Default: csv
        outpath: str
            Path to a wanted output folder. Default None, gets the variable 'outpath'
            previously defined
        outname : str

        withcoord : bool

        coordpath : str, optional

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
    
    def georef(self, var, tag, coordpath, crs = 'epsg:4326', outpath = None,
               fname = None, setindex = False, dropcoord = True, **kwargs):
        """
        Writes a dataframe selected from the ones created in any OGR data
        source supported by Fiona. By default an ESRI shapefile is written.
        
        Parameters
        ----------
        var : str
            As defined in get_df
        tag : str
            As defined in get_df
        coordpath : str
            Path to a .csv file with columns 'row', 'column', self.info['id'], 'x', 'y'
        crs : str, default 'epsg:4326'
            Coordinate reference system. The value can be anything accepted
            by pyproj.CRS.from_user_input()           
        outpath : str, default None
            Path to a wanted output folder. Default: variable 'outpath'
            defined previously
        fname : str, default None
            Name of the file to be written. If None, the tag value will be used
            as name. Here it is needed to specifiy the extension of the file if
            different from '.shp'.
        setindex : bool, default False
            Option to set the dataframe index also in the file created. Default
            is False as it takes lots of time.
        dropcoord : bool, default True
            Drop the coordinates from the DataFrame columns.
        **kwargs : 
            Keyword args to be passed to geopandas.GeoDataFrame.to_file(). For
            example, the file format
        """
        import geopandas as gp
        from shapely.geometry import Point
        
        start = time.time()
        outpath = self.set_outpath(outpath)
        fname = fname if fname else f'{tag}.shp'
        #Get the coordinates from the CSV file and merge them with the dataframe
        coord = pd.read_csv(coordpath)
        coord = self.insert_ind(coord, coord['row'], coord['column'], name = self.info['id'])
        coord = coord.loc[:, (self.info['id'], 'X', 'Y')]
        tool = pd.merge(coord, self.get_df(var, tag), on = self.info['id'])
        #Create the GeoDataFrame
        points = [Point(x,y) for x,y in zip(tool.X,tool.Y)]
        geodf = gp.GeoDataFrame(tool, geometry = points)
        geodf.set_crs(crs, inplace = True)
        if setindex: geodf.set_index(self.info['id'], inplace = True)
        if dropcoord: geodf.drop(['X', 'Y'], axis = 1, inplace = True)
        #Save the GeoDataFrame in the format selected. Default: ESRI Shapefile
        geodf.to_file(f'{outpath}/{fname}', index = setindex, engine = 'fiona', **kwargs)
        end = time.time()
        print(f'Shapefile saved in {outpath} as {fname}')
        print(f'Elapsed time: {round(end-start, 2)} s')
    
    #General functions
    #-----------------
    
    def double_cond(self, c1, c2):
        cond = c1
        for i in range(len(c2)):
            cond = c1 & c2[i]
        return cond
        
    def find_SPcol(self, col, ind = None, indname = False):
        names = [ind] if indname else []
        for name in col:
            if name.find('SP') != -1:
                names += [name]
        return names
    
    def get_df(self, var, tag):
        """
        Returns the selected dataframe from the RechargeCalc object
        
        var (str): 'input', 'recharge'
        tag (str):
         if var: input
         - 'ind', 'irr', 'urb'
         if var: recharge
         - 'rmeteo', 'rirr', 'rurb', 'rtot'
        """
        #Get the df
        if var == 'input':
            df = self.input[tag]
        elif var == 'recharge':
            df = self.recharges[tag]
        return df
        
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
    
    def replicate_columns(self, x, toreplicate):
        print(f"Replicating columns to create a dataframe of {self.info['nSP']} SPs")
        # x = len(self.find_SPcol(irr))
        y = self.info['nSP']
        k, d = divmod(y/x, 1)
        k = round(k)
        tool = toreplicate.set_index('indicatore').copy()
        concatenated = pd.concat([tool[self.find_SPcol(tool)]]*k, axis = 1)
        sps = self.find_SPcol(toreplicate, 'indicatore', True)
        df = toreplicate.loc[:, sps[0:round(x*d)+1]].copy()
        joined = concatenated.join(df.set_index(self.info['id']), on = self.info['id'], rsuffix = '_new')
        cc = [f'SP{i+1}' for i in range(len(joined.columns))]
        joined.columns = cc
        joined = joined.reset_index(level = 0)
        replicated = joined.copy()
        return replicated

        
    def sel_recharge(self, meteo, irr, urb):
        self.conditions['meteo'] = meteo
        self.conditions['irr'] = irr
        self.conditions['urb'] = urb
    
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
    
    def set_outpath(self, outpath = None):
        if not outpath:
            outpath = self.paths['outpath'] if 'outpath' in self.paths else self.paths['input_folder']
        else:
            self.paths['outpath'] = outpath
        return outpath
    
    def save(self, outpath):
         #Save a pickle of the object
         pass