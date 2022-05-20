# -*- coding: utf-8 -*-
"""
Created on Fri May 20 17:35:29 2022

@author: paolo
"""

# %% Setup

import pandas as pd
import numpy as np
import os
import glob
from PIL import Image

os.chdir('C:/E-OBS-SWB2')

# %% Functions 

def get_MMDD(m, d):
    if m < 10:
        m = f'0{m}'
    if d < 10:
        d = f'0{d}'
    return m, d

from Python.custom_functions import save_ArcGRID

# %% Get the files path

folderpath = './Data/CCscenarios/Tesi_MM/ichec-rca4/excel_files'

fls = glob.glob(f'{folderpath}/*.xlsx')

outpath = './Export/ASCII/CCscenarios/stefano_ichec-rca4'

#Definisco il range di date della proiezione
#Per automatizzarlo lo potrei leggere da 
dd = pd.date_range("2006-01-01", "2100-12-31")
idx = np.where(np.invert((dd.month == 2) & (dd.day == 29)))
dd = dd[idx]

#A noi interessa la proiezione dal 2019 al 2100 (tra 2014 e 2018 abbiamo i dati E-OBS),
#ma inserisco 2014 per avere corrispondenza di area su tutta la serie,
#nel caso SWB2 dia problemi utilizzando due aree diverse
wd = pd.date_range("2014-01-01", "2100-12-31")
idx_date = np.where(dd.isin(wd))[0]

columns = ['Precipitation (in)', 'TMAX (F)', 'TMIN (F)']
outfolds = ['precip', 'tmax', 'tmin']
outvars = ['PRCP', 'TMAX', 'TMIN']

lookuptable = pd.read_csv('./Data/CCscenarios/stazioni_CCscenarios_3st.csv', sep = ';')

# %% Raster operations

im = np.array(Image.open('./Data/CCscenarios/geodata/Voronoi_3CCstazioni_5000.tif'))

stcodes = np.unique(im)
for stcode in stcodes:
    im[np.where(im == stcode)] = stcode*100

#Con Voronoi 5000
# xll = 469649.913
# yll = 5013284.782
# cellsize = 5000
# nodata = -9999
#EPSG:32632 per importarlo in SWB2

#Voronoi 5000 WGS84
xll = 8.6103017689999994
yll = 45.2781122179999969
cellsize = 0.1
nodata = -9999
#WGS84 per importarlo su SWB2

# %% ASCII creation

import time
start = time.time()

for fl in fls:
    station = fl.split('_')[-1].split('.')[0]
    df = pd.read_excel(fl)
    df = df.loc[idx_date, columns]
    df.index = dd.date[idx_date]
    df[columns[0]] = df[columns[0]] * 25.4 #mm
    df[columns[1]] = (df[columns[1]] - 32) * 5/9 #°C
    df[columns[2]] = (df[columns[2]] - 32) * 5/9 #°C
    df.columns = [f'{i}_{station}' for i in df.columns]
    if fl == fls[0]:
        tool = df.copy()
    else:
        tool = tool.join(df)

for i, var in enumerate(columns):
    pos = []
    for col in tool.columns:
        if col.find(var) != -1:
            pos += [col]
    tool2 = tool[pos].copy()
    tool2.columns = [col.split('_')[-1] for col in tool2.columns]
    savepath = f"{outpath}/{outfolds[i]}"
    if not os.path.exists(savepath):
        os.makedirs(savepath)
    for day in range(len(tool2.index)):
        grid = im.copy()
        for station in tool2.columns:
            stcode = 100 * lookuptable.loc[np.where(lookuptable['station'] == station), 'stcode'].values[0]
            grid[np.where(grid == stcode)] = tool2.iloc[day,:][station]
            MM, DD = get_MMDD(tool2.index[day].month, tool2.index[day].day)
            fname = f'{savepath}/{outvars[i]}_{tool2.index[day].year}_{MM}_{DD}.asc'
            save_ArcGRID(pd.DataFrame(grid), fname, xll, yll, cellsize, nodata)

end = time.time()
print(f'Creazione dei file ASCII fino al 2100: {round(end-start, 2)} s')
#26 minuti

#grid può essere salvato in una variabile 3d e poi inserito in un netCDF
#salvare poi un netcdf per ogni anno come necessario a SWB2

# %% Add the files from E-OBS

from Python.EOBSobject import EOBSobject

st = time.time()

inpath = './Data/E-OBS'
var = ['rr', 'tx', 'tn']
outpatheobs = f'{outpath}/eobs'
coord = {'lon': [xll, 9.6103017689999994],
          'lat': [yll, 45.6781122180000025]}
coord = pd.DataFrame(coord)
start = 2014
end = 2018
outnames = ['precip', 'tmax', 'tmin']

for i, v in enumerate(var):
    f = EOBSobject(inpath, v, outpatheobs, swb2 = True)
    f.load()
    f.set_outname(outnames[i])
    f.cut_spacetime(coord, start, end, saveformat = 'arcgrid')
    f.close_netcdf()

nd = time.time()

print(f'Creazione dei file ASCII da E-OBS: {round(nd-st, 2)} s')
#41.7 s

# %% zip up everything

import shutil
shutil.make_archive(outpath, 'zip', outpath)
