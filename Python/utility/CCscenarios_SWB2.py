# -*- coding: utf-8 -*-
"""
Climate change scenarios creation for SWB2

- Load the raster with the number of the station as a dataframe
    (be careful on using the same numbers)
- Load the climate change scenarios
- For each climate change scenario:
    For each station:
        - Open the station dataframe
        - Get the ith day
        - Assign the value to the corresponding tile of the TIFF

@author: paolo
"""

import glob
import numpy as np
import os
import pandas as pd
from PIL import Image

#Set the working directory
os.chdir("C:/E-OBS-SWB2")

from Python.custom_functions import save_ArcGRID

def get_MMDD(m, d):
    if m < 10:
        m = f'0{m}'
    if d < 10:
        d = f'0{d}'
    return m, d

#Load the GeoTIFF and set it as a numpy.array
#im = np.array(Image.open('./Data/CCscenarios/geodata/ModelMI_ID_CCstazioni.tif'))

#Fare una prova su qualche giorno per vedere se SWB2 riesce a calcolare con 
# i dati input ritalgiati sull'area del modello, altrimenti ricreare i file
# ASCII partendo dal diagramma di Voronoi creato a monte
im = np.array(Image.open('./Data/CCscenarios/geodata/Voronoi_CCstazioni_10000.tif'))
#ritagliare le righe sopra e sotto il modello, non servono
#poi riottenere le coordinate del punto bottomleft


#Faccio diventare l'array con 336 righe e 660 colonne
# im = np.c_[im, im[:, -1]]
# im = im[:-1, :]
stcodes = np.unique(im)[1:]

lookuptable = pd.read_csv('./Data/CCscenarios/stazioni_CCscenarios.csv', sep = ';')

#Set all the variables needed for save_ArcGRID
#Con ModelMI_ID
xll = 1476150
yll = 5017150
cellsize = 100
nodata = -9999

#Con Voronoi 10000
xll = 469709.677
yll = 4979455.633
cellsize = 10000
nodata = -9999
#EPSG:32632 per importarlo in SWB2

#Set up the folders in which the climate change scenarios are stored
#Example:
# CCsc1 > station1.csv, station2.csv
# CCsc2 > station1.csv, station2.csv

fileformat = 'txt'
filename = 'downscaled_scenario'
numst = 47 #number of stations
outpath = './Export/ASCII/CCscenarios'

start = "2014-01-01"
end = "2100-12-31"

dd = pd.date_range("2006-01-01", "2100-12-31")
idx = np.where(np.invert((dd.month == 2) & (dd.day == 29)))
dd = dd[idx]

wd = pd.date_range("2014-01-01", "2100-12-31")
idx_date = np.where(dd.isin(wd))[0]

#ci servono temperatura massima e minima, al momento invece abbiamo solo la media
variable = ['pav', 'tas']
outfolds = ['prcp', 'tmean']
outvars = ['PRCP', 'TMEAN']

#Con tmax e tmin
# outfolds = ['prcp', 'tmax', 'tmin']
# outvars = ['PRCP', 'TMAX', 'TMIN']
k = 0
for var in variable:
    outfold = outfolds[k]
    outvar = outvars[k]
    folders = [x[0] for x in os.walk(f'./Data/CCscenarios/{var}')]
    files = []
    for folder in folders:
        files += [glob.glob(f'{folder}/*.{fileformat}')]
        tool = []
        for f in files:
            if len(f) != 0:
                tool += [f]
        files = tool
        del tool
        for filepaths in files:
            filepath = next((x for x in filepaths if x.find(filename) != -1),'none')
            filepath = filepath.replace("\\", "/")
            model = filepath.split('/')[-2]
            scenario = filepath.split('/')[-3]
            
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            if not os.path.exists(f"{outpath}/{scenario}"):
                os.makedirs(f"{outpath}/{scenario}")
            if not os.path.exists(f"{outpath}/{scenario}/{model}"):
                os.makedirs(f"{outpath}/{scenario}/{model}")
            if not os.path.exists(f"{outpath}/{scenario}/{model}/{outfold}"):
                os.makedirs(f"{outpath}/{scenario}/{model}/{outfold}")
            savepath = f"{outpath}/{scenario}/{model}/{outfold}"
            
            f = np.loadtxt(filepath, skiprows = 1, usecols = range(1, numst+1))
            columns = open(filepath, 'r').readline().strip().replace('"', '').split(' ')
            
            #estrai solo le righe necessarie (ottieni index tra data inizio e fine e usa quello)
            f = f[idx_date, :]
            #associa righe a colonne
            # df = pd.DataFrame(f, columns = columns)
            #ottieni i codici delle stazioni necessari dal file raster 
            # pd.DataFrame(columns).isin(stcodes)
            #giorno per giorno, prendi una copia del raster e sostituisci il codice
            #della stazione con il valore della stazione in quel giorno
            for n, day in enumerate(f):
                raster = im.copy()
                for stcode in stcodes:
                    #nel caso in cui ci fossero dati esattamente uguali ai numeri inseriti, verrebbero sostituiti
                    #in questo caso non succede perché non sono mai esattamente 9.000 ecc. i dati che usiamo,
                    #ma al fine di migliorare il codice questa issue dovrà essere affrontata
                    st = lookuptable.loc[np.where(lookuptable['stcode'] == stcode)[0], 'columns'].values[0]
                    # st = [lookuptable['columns'][i] for i in range(len(lookuptable['columns'])) if lookuptable['stcode'][i] == stcode][0]
                    st_idx = [i for i in range(len(columns)) if columns[i] == st][0]
                    raster[np.where(raster == stcode)] = day[st_idx]
                    
                    #Salvare un ASCII per ogni giorno vuol dire occupare tantissima memoria
                    #Alternativa 1:
                    #Generare invece un array 3d [giorni, righe, colonne]
                    #Poi salvarlo in un netcdf
                    #Alternativa 2:
                    #Generare uno scenario per volta, poi zipparlo
                    #Alternativa 3:
                    # usare Voronoi 10000, ogni ASCII pesa 800 byte, uno scenario circa 50 byte. Onesto
                    MM, DD = get_MMDD(dd[idx_date][n].month, dd[idx_date][n].month)
                    fname = f'{savepath}/{outvar}_{dd[idx_date][n].year}_{MM}_{DD}.asc'
                    save_ArcGRID(pd.DataFrame(raster), fname, xll, yll, cellsize, nodata)
    k += 1



# for ccpath in cclist:
#     stations = glob.glob(f'{ccpath}/*.{fileformat}')
#     for day in range(start, end): #get the row
#         imarray = np.array(im)
#         for stpath in stations:
#             stcode = 1 #get this from the filename and by interrogating a database of station name+code
#             station = pd.read_csv(stpath)
#             imarray[imarray == stcode] = station.loc[day, 'value'] #insert the value in the raster
#         fname = ''
#         save_ArcGRID(imarray, fname)
