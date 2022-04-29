# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 17:29:21 2022

@author: paolo
"""

from PIL import Image
import numpy as np

# Load .tif file and turn it into a numpy.array
im = Image.open('C:/E-OBS-SWB2/Data/Shp/ModelMI_ID_stazioni.tif')
# im.show()
imarray = np.array(im)
