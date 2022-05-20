# EOBSobject class guide

## Introduction

E-OBS data are bla bla...
Link to the data
Explain NetCDF format
Link to the NetCDF project page

## Code

### Needed modules and classes

### Setup

Import the necessary modules and set up the working directory
```python
import os
import pandas as pd
os.chdir('C:/E-OBS-SWB2')
```

Import the class
```python
from Python.EOBSobject import EOBSobject
```
`Python.` is only needed if the EOBSobject.py file is in another folder as in this repository. If the file is in the same folder as the main only `from EOBSobject` is needed.

### 1. Initialize the class

#### Define the variables
Path to the folder where the E-OBS data are stored: `inpath`. *This is a required parameter*
Path to the folder where to store the results: `outpath`
E-OBS code (lowercase) of the variable needed: `var` *This is a required parameter*
`outpath` can be a path to a custom folder or direct path to the model folder, as the example below.
`outpath` is set as `inpath` if the field is left untouched.
```python
inpath = './Data/E-OBS'
outpath = './Model/swb2_MODELMI/climate_ncfile'
var = 'rr' #daily precipitation sum
```
All the E-OBS variable codes are provided at: https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php
The files in the inpath folder which are not the one you want to load shouldn't have the name of the variable (var) at the start of their file name. The code could identify them as the E-OBS file and thus not work. Example: *rr_otherfile.pdf* shouldn't be in the folder.

#### Create the object
```python
f = EOBSobject(inpath, var, outpath)
```

If you want to provide directly the path to a single file, set `folder` to `False` when creating the object.
```python
inpath = './Data/E-OBS/file.nc'
f = EOBSobject(inpath, var, folder = False)
```

The `swb2` option has to be set to True if the output needs to be used for SWB2.
```python
f = EOBSobject(inpath, var, outpath, swb2 = True)
```

#### Load the netcdf file
```python
f.load()
```

### 2. Generate netcdf files

You can produce netcdf files starting from the E-OBS data.

#### 2.1 Cut in space

Provide the extreme coordinates of the area you want to cut, then call `cut_space()`.
```python
coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)
f.cut_space(coord)
```

It is possible to keep more cells than the area you want to cut, for example 1 cell to 
```python
f.set_fname(f'{outpath}/rr_morecells.nc')
#just to provide a custom name to distinguish the files, not needed for the code to work
f.cut_space(coord, contourcell = 2)
```

#### 2.2 Cut in time

To cut in time, provide the start and end years of the period you choose, then call `cut_time()`.
```python
start = 2014
end = 2018
f.cut_time(start, end)
```

The default option will generate one file for each year. If you want to generate a single file, set `option` as 'bundle'.
```python
f.cut_time(start, end, option = 'bundle')
```

With the `bundle` option, you can also cut between given days. You need to set day as `True` and provide start and end as `datetime.date` objects.
```python
from datetime import date
start_day = date(2011, 12, 30)
end_day = date(2017, 7, 15)
f.cut_time(start_day, end_day, option = 'bundle', day = True)
```

#### 2.3 Cut in space and time

All options available for `cut_space` and `cut_time` are also available for `cut_spacetime`.
```python
f.cut_spacetime(coord, start, end)
```

#### 2.4 Keep the raw file

You can also save the file as it is. This will only change the metadata or other things as the name of the main variable. However, this method is memory consuming and may not work on your laptop as it is.
```python
f.save_netcdf(method = 'raw')
```

### 3. Generate daily ArcGRID files

Specify the method: `'cut_space'`, `'cut_time'`, `'cut_spacetime'`. \
Provide the necessary information: `coord`, `start`, `end`. All the possible information you can provide to obtain a NetCDF you can also apply to this method (for example, the `contourcell` parameter)
```python
f.save_arcgrid('cut_spacetime', coord, start, end)
```
### 4. Perform the operation on multiple E-OBS .nc files

Set up a list of variable tags, then perform a for loop for the different variables.
```python
var = ['rr', 'tn', 'tx']
for v in var:
    f = EOBSobject(inpath, v, outpath)
    f.load()
    f.cut_spacetime(coord, start, end)
    f.close()
```
Outpath and inpath could be provided in lists as well, if you want to refer to single files each time in different folders.