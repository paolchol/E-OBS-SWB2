'''
Modify (i.e., delete the last 4 SP) the following input files from last MODEL-MI version:

- ricarica_irrigua
- extractions
- irrigua_speciale

To be used in RechargeCalc.py
'''
import pandas as pd
import os

#set directory
dir_in = 'C:/Users/user/OneDrive - Politecnico di Milano/GitHub/E-OBS-SWB2/Data/Calcolo_ricarica_totale/'
dir_out = 'C:/Users/user/OneDrive - Politecnico di Milano/GitHub/E-OBS-SWB2/MODEL-MI_2024/Data_modified/'

#open files as dataframes
ext = pd.read_csv(os.path.join(dir_in, 'extractions.csv'))
ric = pd.read_csv(os.path.join(dir_in, 'ricarica_irrigua.csv'))
irr_sp = pd.read_csv(os.path.join(dir_in, 'rirrigua_speciale.csv'))

#Drop last SPs
SPs = ['SP17', 'SP18', 'SP19', 'SP20']
new_ext = ext.drop(SPs, axis=1)
new_ric = ric.drop(SPs, axis=1)
new_irr_sp = irr_sp.drop(SPs, axis=1)

#Save new dfs
new_ext.to_csv(os.path.join(dir_out, 'NEWextractions.csv'))
new_ric.to_csv(os.path.join(dir_out, 'NEWricarica_irrigua.csv'))
new_irr_sp.to_csv(os.path.join(dir_out, 'NEWrirrigua_speciale.csv'))