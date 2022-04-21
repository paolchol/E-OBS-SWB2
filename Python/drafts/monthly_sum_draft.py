# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 17:05:33 2022

@author: paolo
"""


#----------------------------------------------------------
#Analyses section

def monthly_sum(self, coord, start, end,
                loncol = 'lon', latcol = 'lat', contourcell = 0):
    res_cs = self.cut_space(coord, False, True,
                            loncol, latcol, contourcell)
    for year in range(start, end+1):
        res_ct = self.cut_time(start, end, save = False, internal = True)
        df = self.get_var('cut')