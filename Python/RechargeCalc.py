# -*- coding: utf-8 -*-
"""
RechargeCalc class definition
Used to calculate the total recharge
It can be called using:
        from Python.RechargeCalc import RechargeCalc

The working directory has to be set in ./E-OBS-SWB2 for this to work

@author: paolo
"""

import os
import numpy as np
import pandas as pd
import glob
import time

#The directory has to be set in ./E-OBS-SWB2 for this to work
from Python.SWB2output import SWB2output
from Python.custom_functions import repeat_list
# from Python.custom_functions import getkeys

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
        #Store the input files inside the object
        self.input = {
            'ind': pd.read_csv(fls[k[0]])
            }
        if (irr): self.input['irr'] = pd.read_csv(fls[k[1]])
        if (urb): self.input['urb'] = pd.read_csv(fls[k[2]])
        
    def addoutpath(self, outpath):
        self.paths['oupath'] = outpath
    
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
        r = list(map(str, rmeteo['nrow']+1))
        c = list(map(str, rmeteo['ncol']+1))
        newc = []
        for i in range(len(r)):
            newc += [int(f'{r[i]}0{c[i]}')]
        rmeteo.insert(0, self.info['id'], newc)
        
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
            inpath = self.paths['input_folder']
            rmeteo.to_csv(f'{inpath}/rmeteo.csv')

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
        ens = time.time()
        print(f'Elapsed time: {round(end-start, 2)} s')
        if export:
            inpath = self.paths['input_folder']
            rirr.to_csv(f'{inpath}/rirr.csv')

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
            inpath = self.paths['input_folder']
            rurb.to_csv(f'{inpath}/rurb.csv')

    def totalR(self, meteopar, irrpar, urbpar):
        #somma le ricariche
        #scrivere in modo che le componenti da solmmare possano essere scelte
        #autonomamente
        #se gli attributi son o già presenti nella classe, prendere quelli,
        #altrimenti chiamate le funzioni
        print('Total recharge dataframe creation')
        start = time.time()
        
        if 'rirr' not in self.recharges:
            self.irrigationR(irrpar['Is'], irrpar['coeffs'])
        if 'rurb' not in self.recharges:
            self.urbanR(urbpar['coeff_urb'])
        end = time.time()
        print('Elapsed time: {round(end-start, 2) s}')
        
# 	def export(..., data = '', fileformat = '.csv'):
# 		#esportare le ricariche prodotte
# 		#data: urban, irrigation, meteoric, total
# 		#check su data e poi esporta il df voluto nel formato voluto
# 		#format: .csv, .txt
# 	def georef(df, coord):
# 		#esportare un dataframe come shapefile assegnandone le coordinate
