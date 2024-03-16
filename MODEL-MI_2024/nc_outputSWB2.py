"""
Process SWB2 output files
Based on Guide_SWB2output.py by @paolochol

Input files (output from SWB2):
- net_infiltration

Class(es):
- EOBS-SWB2/Python/SWB2output.py
"""

import os
import netCDF4 as nc
import numpy as np
import pandas as pd

#Call the class
os.chdir('C:/Users/HP/Documents/GitHub/E-OBS-SWB2')
from Python.SWB2output import SWB2output

#Path to net infiltration netCDF file produced by SWB2
swb2path = "C:/Users/HP/OneDrive - Politecnico di Milano/SWB2/MODEL-MI/output/R1_net_infiltration__2019-01-01_2022-12-31__338_by_660.nc"

#Create the object
#nobj = nc.Dataset(swb2path)
f = SWB2output(swb2path)
f.metadata #print all the metadata available
print(f.metadata['end_date']) #print one variable of the metadata

#Define stress periods, sum over them and assign sum to a new variable
#Stress period definition
SP1 = 90   #days, 01/01 - 31/03
SP2 = 76   #days, 01/04 - 15/06
SP3 = 92   #days, 16/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]
SPs = np.cumsum(SPs)
SPsum = f.SP_sum(SPs)

f.close()

#Obtain a single dataframe with index, row, column, SP columns
f.obtain_df_SP(index = 'indicatore')
dfSP = f.results['dfSP']
