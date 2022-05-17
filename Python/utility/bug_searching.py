# -*- coding: utf-8 -*-
"""
Prove e test per trovare il bug nella creazione della ricarica totale

@author: paolo
"""

#Prova lanciando RechargeCalc

# %% Setup
import os
import numpy as np
os.chdir('C:/E-OBS-SWB2')
from Python.RechargeCalc_bugfix import RechargeCalc

#Define the variables
startyear = 2014
endyear = 2018
cell_area = 100*100 #m2
bugtrials = "./Data/Calcolo_ricarica_totale/bugtrials"
#Path to the SWB2 output
swb2path = f"{bugtrials}/VersioneFINALE_net_infiltration.nc"
#Path to the input .csv files folder
inputpath = "./Data/Calcolo_ricarica_totale"

r = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore', nSP = 20)
r.load_inputfiles(urb = False)

# %% Calcolo le ricariche

#ottengo rmeteo
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]
rmeteo3d = r.meteoricR(SPs)
rmeteo = r.get_df('recharge', 'rmeteo')

#ottengo rirr
coeffs = {
    'E': 0.3,  
    'R': 0.05, 
    'RISP': 1, 
    'P': 1     
    }
spath = f'{inputpath}/rirrigua_speciale.csv'
r.irrigationR(coeffs, spath)
rirr = r.get_df('recharge', 'rirr')

# %% Carico i file ottenuti da Stefano

import pandas as pd

Smeteo = pd.read_excel(f"{bugtrials}/tabelle_ricaricameteorica.xlsx")
Sirr = pd.read_excel(f"{bugtrials}/tabelle_ricaricairrigua.xlsx")

#sostituisco indicatore
def insert_ind(df, r, c, pos = 0, name = 'none'):
    name = 'indicatore' if name == 'none' else name
    newc = []
    for i in range(len(r)):
        newc += [f'{r[i]}X{c[i]}']
    if (name not in df.columns):
        df.insert(pos, name, newc)
    else:
        df[name] = newc
    return df

Smeteo = insert_ind(Smeteo, Smeteo['row'], Smeteo['column'])
Sirr = insert_ind(Sirr, Sirr['row'], Sirr['column'])

# %% Compare

#Come funziona .isin():
f = pd.DataFrame({'id': [1, 2, 3, 4, 5, 1]})
d = [1, 2]
f.isin(d)

def find_SPcol(col, ind = None, indname = False):
    names = [ind] if indname else []
    for name in col:
        if name.find('SP') != -1:
            names += [name]
    return names

rmeteo.sort_values('indicatore', inplace = True)
rirr.sort_values('indicatore', inplace = True)
Smeteo.sort_values('indicatore', inplace = True)
Sirr.sort_values('indicatore', inplace = True)

#Differenze tra ricariche irrigue
sum(rirr[find_SPcol(rirr)].values == Sirr[find_SPcol(Sirr)].values)
#Le differenze sono nello stesso numero in ogni SP con irrigazione presente

t1 = rirr.loc[Sirr['SP18 [m/sec]'] != rirr['SP18'], ['indicatore', 'distretto', 'SP18']]
t2 = Sirr.loc[Sirr['SP18 [m/sec]'] != rirr['SP18'], ['indicatore', 'NOME_COM', 'SP18 [m/sec]']]
t3 = t1.iloc[:, 2] - t2.iloc[:, 2]

t1 = rirr.loc[Sirr['SP3 [m/sec]'] != rirr['SP3'], ['indicatore', 'distretto', 'SP3']]
t2 = Sirr.loc[Sirr['SP3 [m/sec]'] != rirr['SP3'], ['indicatore', 'NOME_COM', 'SP3 [m/sec]']]
t3 = t1.iloc[:, 2] - t2.iloc[:, 2]
t3.mean()
[t1.iloc[:,2].mean(), t2.iloc[:, 2].mean()]
#Per SP3 e SP18, le differenze partono dall'ordine di 10^(-12), 5 ordini di grandezza inferiori al valore medio
#Non è perciò una differenza che faccia pensare a un errore di codice, 
#al più a differenze nella considerazione delle celle nel calcolo dell'area
#L'errore non è nella ricarica irrigua


locate1 = rmeteo['indicatore'].isin(Smeteo['indicatore'])
locate2 = Smeteo['indicatore'].isin(rmeteo['indicatore'])
sum(rmeteo.loc[locate1, find_SPcol(rmeteo)].values == Smeteo.loc[locate2, find_SPcol(Smeteo)].values)

t1 = rmeteo.loc[locate1, ['indicatore', find_SPcol(rmeteo)[0]]]
t2 = Smeteo.loc[locate2, ['indicatore', 'NOME_COM', find_SPcol(Smeteo)[0]]]
t3 = t1['SP1'] != t2['SP1']

[t1.iloc[:,1].mean(), t2.iloc[:,2].mean()]

rmeteo3d = r.meteoricR(SPs) #return rmeteo3d
#SP1, riga 100, colonna 101
print(rmeteo3d[0, 100, 102])
#SP1, indicatore 100X100: python
print(t1.loc[t1['indicatore'] == '100X101', 'SP1'].values[0])
#SP1, indicatore 100X100: excel
print(t2.loc[t2['indicatore'] == '100X101', 'SP1 [m/s]'].values[0])

#sono tre valori diversi... come è possibile?
#t1: significa che la trasposizione da matrice a lista di righe non ha funzionato
#t2: c'è qualcosa che non va nello script su R

#fare prove su un solo SP, poi dopo generalizzare

#Controllare:
# - definizione dell'indicatore rXc
# - creazione di rmeteo
# - creazione di rirr
# - calcolo di SP_sum


#Confronto con output R

# f = open(f'{bugtrials}/R_ASCII/net_infiltration_2015_SP4.asc', 'r')
# f.read()

Rascii = np.genfromtxt(f'{bugtrials}/R_ASCII/net_infiltration_2015_SP4.asc', dtype = None)
#trasforma in mm
Rascii = Rascii*0.0254/(60*60*24*107)
diff = Rascii - rmeteo3d[7,:,:]

#ci sono delle differenze nell'ordine di grandezza dei dati
#serve controllare che R e python facciano le somme correttamente
#le differenze si trovano nelle celle (numerazione di Python)
[4, 1]
[6, 2]
[26, 14]
[28, 14]
[29, 14]

# %% verifica della correttezza della somma

from Python.SWB2output import SWB2output

swb2 = SWB2output(swb2path)
SPscum = np.cumsum(SPs)
SP_sum_in = swb2.SP_sum(SPscum, units = 'inches', retval = True)
SP_sum_ms = swb2.SP_sum(SPscum, units = 'ms', retval = True)
var = swb2.netCDF['net_infiltration'][:,:,:]


