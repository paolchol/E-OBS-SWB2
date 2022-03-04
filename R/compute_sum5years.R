#Script to compute the sum of the variable
# computed by SWB2 over the whole period considered
#Author: Paolo Colombo

#Works with the variables:
# - 'net_infiltration'
# - 'gross_precipitation'

# Setup -------------------------------------------------------------------

#Inserisci il percorso alla cartella principale
setwd('c:/E-OBS-SWB2')

#Inserisci il percorso alla cartella dove vuoi che le somme vengano salvate
#outpath = "./Export/ASCII/RMeteo_tot" #non mettere l'ultima /
outpath = "./Export/ASCII/Sums"

#Inserisci il percorso alla cartella dove hai salvato gli output di SWB2
#In quella cartella, metti solo i file net_infiltration
#inpath = "./Data/SWB2_output/" #metti l'ultima /
inpath = "./Data/SWB2_output/gross_precipitation/"

#Inserisci il nome della variabile
variable = 'gross_precipitation'

# Library and custom function needed --------------------------------------

library(ncdf4)

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

# Sum over the whole time period ------------------------------------------

fls = list.files(inpath, pattern = '*.nc')

for (i in seq_len(fls)){
  
  f = nc_open(paste0(inpath, fls[i]))
  varsize <- f[["var"]][[variable]][["size"]]

  df = array(0, dim = c(varsize[1], varsize[2]))
  for(j in 1:varsize[3]) df = df + ncvar_get(f, variable, start = c(1, 1, j), count = c(varsize[1], varsize[2], 1))
  df = df*0.0254 #meters
  
  size = round(ncvar_get(f, 'x')[2] - ncvar_get(f, 'x')[1])
  xll = ncvar_get(f, 'x')[1] - size/2
  yll = ncvar_get(f, 'y')[length(ncvar_get(f, 'y'))] - size/2
  
  fname = paste0(outpath, '/', strsplit(fls[i], '_')[[1]][1], "_", variable, "_sum_tot.asc")
  save_ArcGRID(round(df, 6), fname, xll, yll, size, -9999)
  
  nc_close(f)
}
