"""
Get sum of infiltration over stress periods from SWB2 output file
Based on Guide_SWB2output.py by @paolochol

Input files (output from SWB2):
- net_infiltration.nc

Class(es):
- EOBS-SWB2/Python/SWB2output.py

"""

import os
import netCDF4 as nc
import numpy as np
import pandas as pd

#Call the class
os.chdir('c:/Users/user/OneDrive - Politecnico di Milano/GitHub/E-OBS-SWB2')
from Python.SWB2output import SWB2output

#Path to net infiltration netCDF file produced by SWB2
swb2path = "C:/Users/user/OneDrive - Politecnico di Milano/SWB2/MODEL-MI/output/R1_net_infiltration__2019-01-01_2022-12-31__338_by_660.nc"

#Check the file
#nc_file = nc.Dataset(swb2path)

#Create the object
f = SWB2output(swb2path)
f.metadata #print all the metadata available

#Define stress periods, sum over them and assign sum to a new variable
#Stress period definition
SP1 = 90   #days, 01/01 - 31/03
SP2 = 76   #days, 01/04 - 15/06
SP3 = 92   #days, 16/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]
#SPs = np.cumsum(SPs)   #not sure if this is needed

#Assign sum to a new variable and save ASCII files
SPsum = f.SP_sum(SPs, outpath = 'c:/Users/user/OneDrive - Politecnico di Milano/SWB2/DataPreparation/SWB2_output')

f.close()

