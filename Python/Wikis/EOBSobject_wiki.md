# EOBSobject class guide

## Introduction

E-OBS data are bla bla...
Link to the data
Explain NetCDF format
Link to the NetCDF project page

## Code

### Needed modules and classes

The needed modules to use the `EOBSobject` class are: **datetime**, **glob**, **netCDF4**, **numpy**, **os** and **pandas**.

### Setup

Import the necessary modules and set up the working directory.
```python
import os
import pandas as pd
os.chdir('C:/E-OBS-SWB2')
```

Import the class.
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
All the E-OBS variable codes are provided at: https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php. \
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

#### Load E-OBS NetCDF file

Run `load()` to load the actual E-OBS NetCDF file. Once loaded, the original file will be always accessible at `f.netcdf` to perform any operation you want that are available thanks to **netCDF4** package. 
```python
f.load()
```
When you finish doing all the procedures you need on the E-OBS file, remember to close the NetCDF file, since if you keep it open it may consume your laptop RAM.
```python
f.close_netcdf()
```

### 2. Generate NetCDF files

You can produce NetCDF files starting from the E-OBS data. You can make a cut in space, in time and both using `cut_space`, `cut_time` and `cut_spacetime` functions. These functions have their specific requirements explained in the sections below. They all share the possibility to directly save the outcome of the operation performed to the path specified when creating the object. This option is set as default (`save = True`) and will save a NetCDF file as default (`saveformat = 'netcdf'`).

#### 2.1 Cut in space

Provide the extreme coordinates of the area you want to cut, then call `cut_space()`.
```python
coord = {'lon': [8.691, 8.929, 9.524, 9.537],
          'lat': [45.611, 45.308, 45.610, 45.306]}
coord = pd.DataFrame(coord)
f.cut_space(coord)
```

It is possible to keep more cells than the area you want to cut. This could be useful for example for future interpolation to avoid "edge effects".
```python
f.set_fname(f'{outpath}/rr_morecells.nc') #just to provide a custom name to distinguish the files, not needed for the code to work
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

All options available for `cut_space` and `cut_time` (such as contourcells) are also available for `cut_spacetime`.
```python
f.cut_spacetime(coord, start, end)
```

#### 2.4 Keep the raw file

You can also save the file as it is. This will only change the metadata or other things as the name of the main variable. However, this method is memory consuming and may not work on your laptop as it is.
```python
f.save_netcdf(method = 'raw')
```

### 3. Generate daily ArcGRID files

You can generate daily ArcGRID files instead of netCDF files. To do this you can use the same functions used before (`cut_space`, `cut_time` and `cut_spacetime`) specifying the `saveformat` parameter as `'ASCII'`. You have to provide the same information and can apply all the possible parameters used to obtain a NetCDF to this method aside from the `'bundle'` option in `cut_time`.\
This mode will also create a folder in which it will store the ASCII files created. You can set the name of the output folder by using the `set_outname()` function. The output file name will be composed as "{namefile}_YYYY_MM_DD". Namefile will be the same as the namefolder. If `swb2` is set as `True` when create the EOBSobject, `namefile` will be namefolder uppercased. 
```python
f.cut_spacetime(coord, start, end, saveformat = 'arcgrid')
```

### 4. Perform the operation on multiple E-OBS .nc files

If you have multiple files to which you want to perform the operations made available by `EOBSobject`, you can set up a simple for loop. For example, you can set up a list of variable tags, then perform a for loop for the different variables. Outpath and inpath could be provided in lists as well, if you want to refer to single files each time in different folders.
```python
#Basic example
var = ['rr', 'tn', 'tx']
for v in var:
    f = EOBSobject(inpath, v, outpath)
    f.load()
    f.cut_spacetime(coord, start, end)
    f.close()
```

### 5. Examples

#### Working example to generate SWB2 ArcGRID input files

```python
import pandas as pd
inpath = './Data/E-OBS'
var = ['rr', 'tx', 'tn']
outpatheobs = './Export/ASCII/eobs'
coord = {'lon': [8.6103017689999994, 9.6103017689999994],
          'lat': [45.2781122179999969, 45.6781122180000025]}
coord = pd.DataFrame(coord)
start = 2014
end = 2018
outnames = ['precip', 'tmax', 'tmin']

for i, v in enumerate(var):
    f = EOBSobject(inpath, v, outpatheobs, swb2 = True)
    f.load()
    f.set_outname(outnames[i])
    f.cut_spacetime(coord, start, end, saveformat = 'arcgrid')
    f.close_netcdf()

```
