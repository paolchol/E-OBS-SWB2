# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 16:23:57 2022

@author: paolo
"""

def calcarea(ind_df, in_rirr, cell_area):
    in_rirr.insert(len(in_rirr.columns),'area', 0)
    for i, distr in enumerate(in_rirr['distretto'], 0):
        area = sum((ind_df['distretto'] == distr) & (ind_df['zona_agricola'] == 1))*cell_area
        in_rirr.loc[i, 'area'] = area
    return in_rirr

caree = pd.read_csv('./Data/Calcolo_ricarica_totale/copia_calcoloaree.csv')
prova1 = pd.read_csv('./Data/Calcolo_ricarica_totale/ricarica_irrigua.csv')

prova1 = calcarea(caree, prova1, 100*100)
