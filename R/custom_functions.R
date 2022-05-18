#Custom functions
# - General
# - SWB2 netCDF files handling

# General -----------------------------------------------------------------

save_ArcGRID <- function(df, fname, xll, yll, size, nodata){
  datafile <- file(fname, open = 'wt')
  on.exit(close(datafile))
  
  writeLines(paste0('ncols         ', dim(df)[1]), con = datafile)
  writeLines(paste0('nrows         ', dim(df)[2]), con = datafile)
  writeLines(paste0('xllcorner     ', xll), con = datafile)
  writeLines(paste0('yllcorner     ', yll), con = datafile)
  writeLines(paste0('cellsize      ', size), con = datafile)
  writeLines(paste0('NODATA_value  ', nodata), con = datafile)
  
  write(df, datafile, ncolumns = dim(df)[1], sep = ' ')
}

leap <- function(year){
  if(year %% 4 == 0 | year %% 400 == 0) return(366) else return(365)
}

# SWB2 netCDF handling ----------------------------------------------------


SPsum <- function(path, SPs, var,
                  starty, endy,
                  outpath = 'none'){
  #path: path to the netCDF file to be imported
  #SPs: stress periods lengths. Defined as cumulative sums
  library(ncdf4)
  
  print(paste0('Performing the sum of ', var, ' over the stress periods provided'))
  f = nc_open(path)
  varsize = f[["var"]][[var]][["size"]]
  if(outpath != 'none'){
    size = round(ncvar_get(f, 'x')[2] - ncvar_get(f, 'x')[1])
    xll = ncvar_get(f, 'x')[1] - size/2
    yll = ncvar_get(f, 'y')[length(ncvar_get(f, 'y'))] - size/2
  }
  period = seq(starty, endy)
  s = k = 1
  var3d = array(0, dim = c(length(period)*length(SPs), varsize[1], varsize[2]))
  for (y in period){
    year = ncvar_get(f, var, start = c(1, 1, s), count = c(varsize[1], varsize[2], leap(y)))
    base = 1
    for (i in seq_len(length(SPs))){
      SP = ifelse((leap(y) == 366), SPs[i]+1, SPs[i])
      sp = year[, , base:SP]
      base = SP+1
      if((leap(y) == 366) & (i == 1)){
        print(paste0('Be careful, when transforming to m/s, consider that SP1 of year ', y, ' has ', dim(sp)[3], ' days'))
      }
      for(j in seq_len(dim(sp)[3])) var3d[k, , ] = var3d[k, , ] + sp[, , j]
      if (outpath != 'none'){
        fname = paste0(outpath, '/', var, '_', y, '_SP', i, '.asc')
        save_ArcGRID(var3d[k, , ], fname, xll, yll, size, -9999)
      }
      k = k + 1
    }
    s = s + leap(y)
  }
  if (outpath != 'none') print(paste0('The ASCII files are saved in ', outpath))
  nc_close(f)
  return(var3d)
}
