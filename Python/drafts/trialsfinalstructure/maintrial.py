# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:24:02 2022

@author: paolo
"""

#Set the path to the script location (may be useful for the main!)
#all the other folders are then relative to the main position
#in my case, the main.py would be in C:/E-OBS-SWB2

#It doesn't work, for these trials this is the wd
import os
os.chdir('C:/E-OBS-SWB2/Python/drafts/trialsfinalstructure')

# import os
# import sys

# abspath = os.path.abspath(__file__)
# dname = os.path.dirname(abspath)
# os.chdir(dname)

# os.chdir(os.path.dirname(sys.argv[0]))

# sys.executable

# %% Options to call other scripts

#Read the other script
#Same operation as "source" in R
exec(open("xdef.py").read())
x #defined in xdef.py
#__name__ is still main by calling xdef.py like this

#This can be useful to clean up the main code by calling a "module" script and
# running that one

#Import the script as a module
#In that script, set if __name__ == '__main__' to exclude parts of code
#Placing if__name__ == 'test1' makes it run the code below when the script is
# imported

import test1

def service_func():
    print('this is service func')
    
if __name__ == '__main__':
    service_func()
    test1.some_func()
    #__name__ printed by some_func() is test1

#Import the script as a module in another folder
os.chdir('C:/E-OBS-SWB2')

import importlib.util
spec = importlib.util.spec_from_file_location("xdef", "./Python/drafts/trialsfinalstructure/xdef.py")
xdef = importlib.util.module_from_spec(spec)
spec.loader.exec_module(xdef)

xdef.xfunc()

#Import the script as a module in a subfolder of the directory
os.chdir('C:/E-OBS-SWB2')
import Python.custom_functions as cf

cf.leap(2017)



