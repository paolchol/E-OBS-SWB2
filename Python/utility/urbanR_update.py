# -*- coding: utf-8 -*-
"""
Trial for the new urbanR

Features needed:

Calculate the recharge due to losses from the extraction pumps
- apply a coefficient to extractions (provided in the extractions.csv input file)

Assign the recharge obtained to cells with specified conditions
- let the user choose the condition to apply
- assign the recharge to the cells that respect the conditions


@author: paolo
"""


import os
import numpy as np

os.chdir('C:/E-OBS-SWB2')

# %% Call the class

from Python.RechargeCalc import RechargeCalc

# %% Class usage

#0. Set up your input files

#1. Initialize the class

#Define the variables
startyear = 2014
endyear = 2018
cell_area = 100*100 #m2
#Path to the SWB2 output
swb2path = "./Data/SWB2_output/VersioneFINALE_net_infiltration.nc"
#Path to the input .csv files folder
inputpath = "./Data/Calcolo_ricarica_totale"

r = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore', nSP = 20)

#Load the input files needed
r.load_inputfiles()

cond = r.urbanR(0.15, False, True,
                col = ['land_cover', 'land_cover', 'zona_urbana'],
                valcol = [123, 124, 1],
                option = [0, 1])

c1 = r.input['ind']['zona_urbana'] == 1
c2 = r.input['ind']['land_cover'] == 123
c3 = r.input['ind']['land_cover'] == 124
sum(c1 & (c2 | c3))
sum(cond)

rurb = r.urbanR(0.15, False, True,
                col = ['land_cover', 'land_cover', 'zona_urbana'],
                valcol = [123, 124, 1],
                option = [0, 1])
