
'''
To just check the netcdf files and change the name of the variables to match SWB2 requirements

'''
from netCDF4 import Dataset
import netCDF4 as nc
import os
import numpy as np

#Path to files
files_dir = 'C:/Users/HP/OneDrive - Politecnico di Milano/SWB2/'
#inpath = os.path.join(files_dir,'DataPreparation/EOBS_object/OutputData/')
inpath2 = os.path.join(files_dir,'MODEL-MI/climate_ncfile/')

#Load and check coords of nc
# prcp_19 = nc.Dataset(os.path.join(inpath,'rr_EOBS_cut_spacetime_2019.nc'), 'r+')
# x_coords = prcp_19.variables['x'][:]
# y_coords = prcp_19.variables['y'][:]
# prcp_19.close()

#Rename variables
# prcp = ['rr_EOBS_cut_spacetime_2020.nc','rr_EOBS_cut_spacetime_2021.nc','rr_EOBS_cut_spacetime_2022.nc']
# tmin = ['tn_EOBS_cut_spacetime_2019.nc','tn_EOBS_cut_spacetime_2020.nc','tn_EOBS_cut_spacetime_2021.nc','tn_EOBS_cut_spacetime_2022.nc']
# tmax = ['tx_EOBS_cut_spacetime_2019.nc','tx_EOBS_cut_spacetime_2020.nc','tx_EOBS_cut_spacetime_2021.nc','tx_EOBS_cut_spacetime_2022.nc']

prcp = ['prcp_EOBS_2019.nc','prcp_EOBS_2020.nc','prcp_EOBS_2021.nc','prcp_EOBS_2022.nc']
tmin = ['tmin_EOBS_2019.nc','tmin_EOBS_2020.nc','tmin_EOBS_2021.nc','tmin_EOBS_2022.nc']
tmax = ['tmax_EOBS_2019.nc','tmax_EOBS_2020.nc','tmax_EOBS_2021.nc','tmax_EOBS_2022.nc']

for p in prcp:
    f = nc.Dataset(os.path.join(inpath2,p), 'r+')
    f.renameVariable(u'rr',u'prcp')
    f.close()
    
for n in tmin:
    f = nc.Dataset(os.path.join(inpath2,n), 'r+')
    f.renameVariable(u'tn',u'tmin')
    f.close()

for m in tmax:
    f = nc.Dataset(os.path.join(inpath2,m), 'r+')
    f.renameVariable(u'tx',u'tmax')
    f.close()