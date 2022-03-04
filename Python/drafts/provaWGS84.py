# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 14:44:17 2022

@author: paolo
"""

outpath = "./Export/ASCII/RMeteo_tot" #not include /
inpath = "./Model/swb2_MODELMI/output/ModelMI_gross_precipitation__2014-01-01_2018-12-31__338_by_660.nc" #include /
# inpath = "./Model/swb2_MODELMI/output/ModelMI_net_infiltration__2014-01-01_2018-12-31__338_by_660.nc" #include /

#fl = glob.glob(f'{inpath}*.nc')

variable = 'gross_precipitation'
# variable = 'net_infiltration'

f = nc.Dataset(inpath)
#net_infiltration = np.ma.getdata(f['net_infiltration'][:,:,:])

df = np.sum(np.ma.getdata(f[variable][:,:,:]), axis = 0)*0.0254 #meters
df = pd.DataFrame(df)

size = round(np.ma.getdata(f['x'][1]).item() - np.ma.getdata(f['x'][0]).item()) #controlla che sia 100
xll = round(np.ma.getdata(f['x'][0]).item()) - size/2
yll = round(np.ma.getdata(f['y'][-1]).item()) - size/2

fname = f"{outpath}/pWGS84_1_{variable}_sum_tot.asc" #To obtain a CSV just change to .csv and run this and the below line
save_ArcGRID(df, fname, xll, yll, size, -9999)

f.close()

# %% Confronto tra gross_precipitation e dati E-OBS originali (pWGS84)

f = nc.Dataset(inpath)
var = np.ma.getdata(f[variable][1:10, :, :])*25.4 #mm
eobs = nc.Dataset(r".\Model\swb2_MODELMI\climate_ncfile\prcp_EOBS_2014.nc")
obsvar = np.ma.getdata(eobs['prcp'][1:10, :, :]) #mm

eobs.close()
