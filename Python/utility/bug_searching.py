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
r.meteoricR(SPs)
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

Smeteo = r.insert_ind(Smeteo, Smeteo['row'], Smeteo['column'])
Sirr = r.insert_ind(Sirr, Sirr['row'], Sirr['column'])

# %% Compare

#Come funziona .isin():
f = pd.DataFrame({'id': [1, 2, 3, 4, 5, 1]})
d = [1, 2]
f.isin(d)

rmeteo.sort_values('indicatore', inplace = True)
rirr.sort_values('indicatore', inplace = True)
Smeteo.sort_values('indicatore', inplace = True)
Sirr.sort_values('indicatore', inplace = True)

#Differenze tra ricariche irrigue
sum(rirr[r.find_SPcol(rirr)].values == Sirr[r.find_SPcol(Sirr)].values)
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

#SP8
Rascii = np.genfromtxt(f'{bugtrials}/R_ASCII/net_infiltration_2015_SP4.asc', dtype = None)
#trasforma in mm
Rascii = Rascii*0.0254/(60*60*24*107)
diff = Rascii - rmeteo3d[7,:,:]

#ci sono delle differenze nell'ordine di grandezza dei dati
#serve controllare che R e python facciano le somme correttamente
#le differenze si trovano nelle celle (numerazione di Python)
check_val = [[4, 1],
             [6, 2],
             [26, 14],
             [28, 14],
             [29, 14]]

# %% Verifica della correttezza della somma

from Python.SWB2output import SWB2output

swb2 = SWB2output(swb2path)
SPscum = np.cumsum(SPs)
SP_sum_in = swb2.SP_sum(SPscum, units = 'inches', retval = True)
SP_sum_ms = swb2.SP_sum(SPscum, units = 'ms', retval = True)
var = swb2.netCDF['net_infiltration'][:,:,:]

#Stress period 8
SP8_start = SP1 + SP2 + SP3 + SP4 + SP1 + SP2 + SP3
SP8_end = SP8_start + SP4

SP8 = var[SP8_start:SP8_end, 4, 1]
sumsp8 = np.sum(SP8, axis = 0)

print(Rascii[4,1] - sumsp8)
print(SP_sum_in[7, 4, 1] - sumsp8)

#C'è un errore nel calcolo che viene effettuato su R

for pos in check_val:
    SP8 = var[SP8_start:SP8_end, pos[0], pos[1]]
    sumsp8 = np.sum(SP8, axis = 0)
    print(f'Difference R: {Rascii[pos[0], pos[1]] - sumsp8}')
    print(f'Difference Python: {SP_sum_in[7, pos[0], pos[1]] - sumsp8}')

#Sistemare il calcolo su R e rifare prova

Rascii_fix = np.genfromtxt(f'{bugtrials}/R_ASCII/fix/net_infiltration_2015_SP4.asc', dtype = None)
for pos in check_val:
    SP8 = var[SP8_start:SP8_end, pos[0], pos[1]]
    sumsp8 = np.sum(SP8, axis = 0)
    print(f'Difference R: {Rascii_fix[pos[0], pos[1]] - sumsp8}')
    print(f'Difference Python: {SP_sum_in[7, pos[0], pos[1]] - sumsp8}')

#Ora i risultati sono uguali
#Provo su un altro stress period
#SP13
Rascii_fix = np.genfromtxt(f'{bugtrials}/R_ASCII/fix/net_infiltration_2017_SP1.asc', dtype = None)
diff = Rascii_fix - SP_sum_in[12, :, :]
sum(sum(diff > 10E-8)) #35453 celle hanno una differenza maggiore di 10^-8
sum(sum(diff > 10E-7)) #0 celle hanno una differenza maggiore di 10^-7

#La somma operata sia da SWB2output.SP_sum sia dalla funzione SP_sum su R è ora corretta

diff = Rascii - SP_sum_ms[7, :,:]


# %% Confronto tra rmeteo
#Generare un rmeteo tramite ASCII su QGIS (primi 2 SP)
#Generare rmeteo tramite meteoricR con RechargeCalc

#load R_gen_rmeteo_inches.shp
#inserisci indicatore con la x
import geopandas as gpd
R_gen = gpd.read_file(f'{bugtrials}/R_gen_rmeteo_inches.shp')
R_gen = r.insert_ind(R_gen, R_gen['row'], R_gen['column'])

#genera rmeteo per i primi due SP
#rmeteo.insert

#rmeteo.join
rmeteoj, rmeteoi = r.meteoricR(SPs, units = 'inches')
rmeteoi.sort_values('indicatore', inplace = True)
rmeteoj.sort_values('indicatore', inplace = True)
for sp in r.find_SPcol(rmeteoi):
    diff = rmeteoi[sp] - rmeteoj[sp]
    print(sp)
    print(f'Numero di celle con diff > 10^-8: {sum(diff > 10E-8)}')
    print(f'Numero di celle con diff > 10^-7: {sum(diff > 10E-7)}')
#I due metodi forniscono lo stesso risultato
#Mantenere quello più efficiente (in teoria è join da quello che dice la documentazione di pandas)

#Confronto con R
#ordina per indicatore
#fai differenza
#conta le differenze maggiori di 10E-8
joined = rmeteoj.join(R_gen.loc[:, ['indicatore','SP1_1', 'SP2_1']].set_index('indicatore'), on = 'indicatore')
filtered = joined[~joined['SP1_1'].isnull()].copy()
filtered['diff1'] = filtered['SP1_1'] - filtered['SP1']
filtered['diff2'] = filtered['SP2_1'] - filtered['SP2']
sum(filtered['diff1'] > 10E-5)
sum(filtered['diff2'] > 10E-5)
#circa 60000 celle sono diverse

# %% Cambiamento nell'assegnamento dell'indicatore 
"""
In meteoricR, quando traformo da griglia a righe/colonne, alla prima colonna
della griglia va assegnato il numero 4, corrispondente alla prima colonna
nello shapefile di punti utilizzato per indicatori.csv

Sommare quindi 4 invece di 1 alle colonne, quando viene inserita la colonna indicatore

Mettere opzione in insert_ind e anche in meteoricR
"""

#Aggiornamento di meteoricR

rmeteo = r.meteoricR(SPs, 'inches', 1, 4, ret = True)
joined = rmeteo.join(R_gen.loc[:, ['indicatore','SP1_1', 'SP2_1']].set_index('indicatore'), on = 'indicatore')
filtered = joined[~joined['SP1_1'].isnull()].copy()
filtered['diff1'] = filtered['SP1_1'] - filtered['SP1']
filtered['diff2'] = filtered['SP2_1'] - filtered['SP2']
print(f"Numero di celle con diff > 10^-7 (SP1): {sum(filtered['diff1'] > 10E-7)}")
print(f"Numero di celle con diff > 10^-7 (SP2): {sum(filtered['diff2'] > 10E-7)}")

print(filtered.loc[filtered['diff1'] > 10E-7, 'diff1'])
#differenze invisibili

# %% Ottimizzazione creazione meteoricR

f = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore', nSP = 153)
f.load_inputfiles(urb = False)
rmeteo160sp = f.meteoricR(SPs, 'inches', 1, 4, ret = True)


import pandas as pd

x = 20
y = 60

k, d = divmod(y/x, 1)
k = round(k)

# sps = r.find_SPcol(rmeteo, 'indicatore', True)

# cc = pd.concat([rmeteo[r.find_SPcol(rmeteo)]]*k, axis = 1).columns

prova = rmeteo.set_index('indicatore').copy()
sps = r.find_SPcol(rmeteo, 'indicatore', True)

concatenated = pd.concat([prova[r.find_SPcol(prova)]]*k, axis = 1)

df = rmeteo.loc[:, sps[0:round(x*d)+1]].copy()
joined = concatenated.join(df.set_index('indicatore'), on = 'indicatore', rsuffix = '_new')
cc = joined.columns
cc = [f'SP{i+1}' for i in range(len(cc))]
joined.columns = cc
joined = joined.reset_index(level = 0)
print(joined.columns)


# %% Confronto tra rtot


#carica tabella rmeteo stefano
#settare come rmeteo dentro rechargecalc

#calcolare rirr
#sommare

#fare confronto tra rtot stefano ed rtot calcolata da RechargeCalc

r.meteoricR(SPs, 'ms', 1, 4)
r.irrigationR(coeffs, spath)
r.totalR()
rtot = r.get_df('recharge', 'rtot')

import pandas as pd
Stot = pd.read_excel(f"{bugtrials}/tabelle_ricaricatotale.xlsx")
Stot = r.insert_ind(Stot, Stot['row'], Stot['column'])
col = Stot.columns
col2 = [c for c in col]
col2[10:] = [f'SP{i+1}' for i in range(20)]
Stot.columns = col2

compare = rtot.join(Stot.set_index('indicatore'), on = 'indicatore', rsuffix = 'S')

col1 = [f'SP{i+1}' for i in range(20)]
col2 = [f'SP{i+1}S' for i in range(20)]

diff = compare['SP2'] - compare['SP2S']

rmeteo = r.get_df('recharge', 'rmeteo')
def op (Stot):
    Stot = r.insert_ind(Stot, Stot['row'], Stot['column'])
    col = Stot.columns
    col2 = [c for c in col]
    col2[10:] = [f'SP{i+1}' for i in range(20)]
    Stot.columns = col2
    return Stot

Smeteo = op(Smeteo)

compare = rmeteo.join(Smeteo.set_index('indicatore'), on = 'indicatore', rsuffix = 'S')

sum(diff > 10E-16)
