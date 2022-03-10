# -*- coding: utf-8 -*-
"""
Notes on how to get information from a netCDF file

@author: paolo
"""

# %% How to get information from a netCDF file

#Variables is a dictionary
#Get the names of the variables
var = f.netCDF.variables #prints all the variables and their metadata
keys = f.netCDF.variables.keys()
print(keys)

def getkeys(dict):
    #got from: https://www.geeksforgeeks.org/python-get-dictionary-keys-as-a-list/
    # return list(dict.keys()) other method
    return [*dict]

keys = getkeys(f.netCDF.variables)
getkeys(f.netCDF.variables)[3] #only the main variable

#Start and end dates

f.netCDF.variables['time'].units #starting date
f.netCDF.variables['time'].shape[0] #time length

def getdates(string, n):
    from datetime import date, timedelta
    datestr = string.split(' ')[2]
    y, m, d = datestr.split('-')
    start = date(int(y), int(m), int(d))
    end = start + timedelta(days = n-1)
    return start, end

string = f.netCDF.variables['time'].units
n = f.netCDF.variables['time'].shape[0]
start, end = getdates(string, n)

#Units
f.netCDF[getkeys(f.netCDF.variables)[3]].units


