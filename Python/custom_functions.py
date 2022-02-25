"""
Collection of custom functions produced in the development of the E-OBS-SWB2 project
Functions are divided in four sections:
0. Generic functions, useful in many contexts
1. E-OBS data handling
2. ArcASCII GRID creation (general)
3. netCDF input creation (SWB2 specific)

author: paolo
"""

# %% 0. Generic functions

def leap(y):
    #input: year (int)
    #output: number of days (int)
    if((y%4 == 0) | (y%400 == 0)):
        return 366
    else:
        return 365

# %% 1. E-OBS data handling

def date_toeobs(y,m,d):
    #Returns how many days have passed from 1950-01-01 (starting counting day for E-OBS data)
    from datetime import date
    start = date(1950,1,1)
    d = date(y,m,d)
    count = d - start #how many days have passed since start
    return count.days

def eobs_todate(x, number = False):
    #Returns the year, month and day corresponding to the number given
    #If x is a plain number (not a variable), "number" must be set to True
    from datetime import date, timedelta
    start = date(1950,1,1)
    if (number): end = start + timedelta(days = x.item())
    else: end = start + timedelta(days = x)
    return end.year, end.strftime('%m'), end.strftime('%d')

def transf(lat_t, lon_t, zoneN, zoneL, var = 'xy'):
    #Transforms the lat/lon data provided in x/y projected UTM coordinates
    import utm
    x = []
    y = []
    for i in lat_t:
        for j in lon_t:
                #The syntax is utm.from_latlon(LATITUDE, LONGITUDE)
                #The return has the form (EASTING, NORTHING, ZONE_NUMBER, ZONE_LETTER)
            xj, yi, _, _ = utm.from_latlon(i, j,
                            force_zone_number = zoneN,
                            force_zone_letter = zoneL)
            x += [xj]
        y += [yi]
    x = x[0:len(lon_t)]
    if(var == 'xy'):
        return x, y
    elif(var == 'x'):
        return x
    elif(var == 'y'):
        return y

# %% 2. ArcASCII GRID creation (general)

def save_ArcGRID(df, fname, xll, yll, size, nodata):
    #df has to be a Pandas DataFrame
    #xll: x coordinate of the left bottom corner (lon)
    #yll: y coordinate of the left bottom corner (lat)
    #size: cell size (m)
    #nodata: value assigned to nodata
    def line_prepender(filename, line):
        with open(filename, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line + '\n' + content)
            #line.rstrip('\r\n') if you want ot remove something from line
    df.to_csv(fname, sep = ' ', header = False, index = False)
    header = f'ncols         {len(df.columns)}\nnrows         {len(df.index)}\nxllcorner     {xll}\nyllcorner     {yll}\ncellsize      {size}\nNODATA_value  {nodata}'
    line_prepender(fname, header)


# %% 3. netCDF input creation (SWB2 specific)

def eobs_todaymet(y):
    #redefine: starting from 1980, adding 0.5
    from datetime import date
    dstart = date(1980, 1, 1)
    estart = date(1950, 1, 1)
    k = dstart - estart
    y1 = y - k.days + 0.5
    return y1