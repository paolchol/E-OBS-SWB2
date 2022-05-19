# RechargeCalc class guide

## Introduction

### What this class does
It calculates a total groundwater recharge starting from three possible components: meteoric, irrigation and urban.
The meteoric recharge has to be provided as a NetCDF output file of SWB2 [link to SWB2 reference].

### 0. Set up your input files

In the folder *./Data/Calcolo_ricarica_irrigua/template* in this repository you can find the templates needed to be able to use the RechargeCalc class. The files are in the CSV format and should remain in that format. To modify them you can use a spreadsheet software as MS Excel or LibreOffice Calc, then save them as CSV being sure that the separator is the comma. Below the instructions on how to fill them. It is important that the original name remains in the final file you create. For example *"indicatori_v1.csv"* is ok, but *"ind_v1.csv"* will not work. The folder where you save the input files should not have files with duplicated names, as the class will then load only the first one (e.g. *"indicatori_v1.csv"* and *"indicatori_v2.csv"* should not be in the same folder).

*indicatori.csv*
| Variable | Description |
| -------- | ----------- |
| FID      | Number of the point. Not used, you can place arbitrary numbers |
| row      | Row of the point. `int` |
| column   | Column of the point. `int` |
| nome_com | Name of the municipality the point falls in. `string`|
| sig_pro  | Abbreviation of the province the point falls in. `string`|
| land_cover | Land cover classification. Not used, you can place arbitrary numbers |
| zona_agricola | 1 if the point is in an agricoltural zone, 0 otherwise. `int` |
| zona_urbana | 1 if the point is in an urban zone, 0 otherwise. `int` |
| distretto | Name of the district the point falls in. `string` |

*ricarica_irrigua.csv*
| Variable | Description |
| -------- | ----------- |
| distretto |  |
| code |  |
| SPi |  |

*rirrigua_speciale.csv*
| Variable | Description |
| -------- | ----------- |
| distretto |  |
| SPi |  |

*extractions.csv*
| Variable | Description |
| -------- | ----------- |
| nome_com |  |
| SPi |  |

*coord.csv*
| Variable | Description |
| -------- | ----------- |
| row |  |
| column |  |
| indicatore |  |
| X |  |
| Y |  |

## Code

### Needed modules and classes

The needed modules to use the `RechargeCalc` class are: glob, numpy, pandas, sys and time. To run `RechargeCalc.georef()` also geopandas is needed. 
The `SWB2output` class is also needed (available in this GitHub repository), which in turn makes it needed to import netCDF4 and a series of functions from *custom_functions.py*. By downloading the repository */Python* folder as it is, the "concatenation" of classes and functions will work automatically by setting the working directory in the parent folder of where you will store the */Python* folder.

### Setup

Import the necessary modules and set up the working directory
```python
import os
import numpy as np
os.chdir('C:/E-OBS-SWB2')
```

Import the class
```python
from Python.RechargeCalc import RechargeCalc
```

### 1. Initialize the class

Define the variables
```python
startyear = 2014
endyear = 2018
cell_area = 100*100 #m2
```

Define the paths
```python
#Path to the SWB2 output
swb2path = "./Data/SWB2_output/1Speranza_netinfiltration.nc"
#Path to the input .csv files folder
inputpath = "./Data/Calcolo_ricarica_totale"
```

Create the object
```python
r = RechargeCalc(swb2path, inputpath, startyear,
                 endyear, cell_area, uniqueid = 'indicatore', nSP = 20)
```

Load the input files needed
```python
r.load_inputfiles()
```

Here you have to select the recharge components you want to consider. If you don't want to consider a component, set its parameter (`meteo`, `urb` or `irr`) as `False`.
```python
r.load_inputfiles(urb = False) #Not consider the urban recharge
```

Once you loaded the input files, you can access them via `r.get_df()` by specifying `'input'` and the desired dataframe tag (`'ind'`, `'irr'` or `'urb'`). The meteoric input file can not be accessed in this way.
```python
r.get_df('input', 'ind') #Extract "indicatori.csv"
```

You can then perform any operation you would on pandas DataFrames for example correct a wrong value provided
```python
r.input['ind'].loc[r.input['ind']['distretto'] == 'Muzza', 'distretto'] = 'MUZZA'
```

### 2. Create meteoric recharge dataframe

Define the stress periods lengths
```python
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = [SP1, SP2, SP3, SP4]
```

Launch `meteoricR()` by providing the SPs, the units wanted and the starting row and column of *"indicatori.csv"* file to correctly associate the cells from the SWB2 output to the cells of the desired area.
```python
r.meteoricR(SPs, 'ms', 1, 4)
```

### 3. Create irrigation recharge dataframe

Define the coefficients needed
coeffs = {
    'E': 0.3,  #Irrigation technique efficiency
    'R': 0.05, #Residual runoff
    'RISP': 1, #1 - fraction of water saved by a change of irrigation technique
    'P': 1     #Percentage of the cell covered by the irrigation
    }
#Path to the input file related to the "special" irrigation district
#Leave it as 'none' if you don't have one
spath = f'{inputpath}/rirrigua_speciale.csv'

r.irrigationR(coeffs, spath)

### 4. Create urban recharge dataframe

coeff_urb = 0.15

r.urbanR(coeff_urb)

### 5. Create total recharge dataframe

#### 5.1 After computing the needed partial recharges
Launch after computing the needed partial recharges. You don't need to have computed all the recharges, just the ones you need.
```python
r.totalR()
```

#### 5.2 Launch directly
Needs the parameters defined before, inside a dictionary or a dataframe as below. The recharges that will be computed are the ones that have its parameter in `load_inputfiles()` set as `True`.
```python
meteopar = {
    'SPs': SPs
    }
irrpar = {
    'coeffs': coeffs,
    'spath': spath
    }
urbpar = {
    'coeff_urb': coeff_urb
    }

r.totalR(meteopar, irrpar, urbpar)
```

### 6. Export

#as .csv
outpath = "./Data/Calcolo_ricarica_totale"
r.export('recharge', 'rtot', outpath = outpath)
#as .shp
#proj: crs in the PROJ4 format. In this case, Monte Mario EPSG:3003
r.georef('recharge', 'rtot', "./Data/Calcolo_ricarica_totale/coord.csv",
         proj = '+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +units=m +no_defs ',
         outpath = outpath)
#If you don't want to keep the X and Y columns in the exported dataframe,
#set dropcoord as True
r.georef('recharge', 'rtot', "./Data/Calcolo_ricarica_totale/coord.csv",
         proj = '+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=1500000 +y_0=0 +ellps=intl +towgs84=-104.1,-49.1,-9.9,0.971,-2.917,0.714,-11.68 +units=m +no_defs ',
         outpath = outpath, outname = 'rtot_noXY', dropcoord = True)


