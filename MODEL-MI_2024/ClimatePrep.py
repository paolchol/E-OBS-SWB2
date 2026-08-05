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
os.chdir('C:/Users/user/Documents/GitHub/E-OBS-SWB2')
from Python.EOBSobject import EOBSobject

# Paths to the input and output data folders
files_dir = 'E:/Neuchatel/Estensione_past/EOBS_object/'
inpath = os.path.join(files_dir,'InputData')
outpath = os.path.join(files_dir,'OutputData')

# %%
# Additional: for a long time series
import xarray as xr
import glob

tx_files = sorted(glob.glob(os.path.join(inpath, 'tx_*.nc')))
tn_files = sorted(glob.glob(os.path.join(inpath, 'tn_*.nc')))
rr_files = sorted(glob.glob(os.path.join(inpath, 'rr_*.nc')))

tx_ds = xr.open_mfdataset(tx_files, combine='by_coords')
tn_ds = xr.open_mfdataset(tn_files, combine='by_coords')
rr_ds = xr.open_mfdataset(rr_files, combine='by_coords')

tx_ds.to_netcdf(os.path.join(inpath, 'tx_1950-2023.nc'))
tn_ds.to_netcdf(os.path.join(inpath, 'tn_1950-2023.nc'))
rr_ds.to_netcdf(os.path.join(inpath, 'rr_1950-2023.nc'))

tx_ds.close()
tn_ds.close()  
rr_ds.close()

# make sure that these are the only datasets in the input folder before running the next cell.

# %%
# Extreme coords of MODEL-MI area
coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)

#Period for aggiornamento MODEL-MI
start = 1953
end = 2022

#Cut E-OBS data in space and time
var = ['rr', 'tn', 'tx']




# %%
for v in var:
    f = EOBSobject(v, inpath, outpath, folder = True, swb2 = True)
    f.load()
    f.cut_spacetime(coord, start, end, contourcell=2)

# %%
