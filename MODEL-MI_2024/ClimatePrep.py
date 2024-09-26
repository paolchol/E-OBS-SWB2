"""
Prepare climate data for area and years of interest --> aggiornamento MODEL-MI 2019-2023
Based on Guide_EOBSoject.py by @paolochol

Input files (E-OBS data from Copernicus):
- daily maximum temperature (TX)
- daily minimum temperature (TN)
- daily precipitation sum (RR)

Class(es):
- EOBS-SWB2/Python/EOBSobject.py

"""

# %% Setup

import os
import pandas as pd
from datetime import date

# %% Call the class
os.chdir('C:/Users/HP/Documents/GitHub/E-OBS-SWB2')
from Python.EOBSobject import EOBSobject

# Paths to the input and output data folders
files_dir = 'C:/Users/HP/OneDrive - Politecnico di Milano/SWB2/DataPreparation/EOBS_object/'
inpath = os.path.join(files_dir,'InputData')
outpath = os.path.join(files_dir,'OutputData')

# Extreme coords of MODEL-MI area
coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)

#Period for aggiornamento MODEL-MI
start = 2019
end = 2022

#Cut E-OBS data in space and time
var = ['rr', 'tn', 'tx']

for v in var:
    f = EOBSobject(v, inpath, outpath, folder = True, swb2 = True)
    f.load()
    f.cut_spacetime(coord, start, end, contourcell=2)