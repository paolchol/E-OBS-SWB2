# -*- coding: utf-8 -*-
"""
Created on Tue May 31 16:10:38 2022

@author: paolo
"""

import numpy as np
import pandas as pd
import os

os.chdir('C:/E-OBS-SWB2')
from Python.EOBSobject import EOBSobject

outpath = './Stefano/eobs_sumSP/'

inpath = './Data/E-OBS'
var = 'rr'
coord = {'lon': [8.6103017689999994, 9.6103017689999994],
          'lat': [45.2781122179999969, 45.6781122180000025]}
coord = pd.DataFrame(coord)
start = 2014
end = 2018

SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]

f = EOBSobject(inpath, var, outpath, swb2 = True)
f.load()
f.SP_sum(SPs, outpath = outpath,
                 coord = coord, start = start, end = end,
                 contourcell = 2,
                 method = 'cut_spacetime',
                 units = 'mm')
f.close_netcdf()

