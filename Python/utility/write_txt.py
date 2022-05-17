# -*- coding: utf-8 -*-
"""
Created on Fri May 13 15:30:31 2022

@author: paolo
"""

h = open("columns.txt", "w")
for col in columns: h.write(col + "\n")
h.close()
