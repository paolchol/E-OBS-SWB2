#Script to compute the sum of the net infiltration
# computed by SWB2 over the whole period considered
#Author: Paolo Colombo

# Setup -------------------------------------------------------------------

#Inserisci il percorso alla cartella principale
setwd('c:/E-OBS-SWB2')

#Inserisci il percorso alla cartella dove vuoi che le somme vengano salvate
outpath = "./Export/ASCII/RMeteo_tot" #non mettere l'ultima /

#Inserisci il percorso alla cartella dove hai salvato gli output di SWB2
#In quella cartella, metti solo i file net_infiltration
inpath = "./Data/SWB2_output/" #metti l'ultima /

# Library and custom function needed --------------------------------------

library(ncdf4)

save_ArcGRID <- function(df, fname, xll, yll, size, nodata){
  datafile <- file(fname, open = 'wt')
  on.exit(close(datafile))
  
  header <- paste0(
            'ncols         ', dim(df)[2],
            '\nnrows         ', dim(df)[1],
            '\nxllcorner     ', xll,
            '\nyllcorner     ', yll,
            '\ncellsize      ', size,
            '\nNODATA_value  ', nodata)
  writeLines(header, con = datafile)
  write.table(df, datafile, sep = ' ', row.names = FALSE,
              col.names = FALSE, quote = FALSE)
}

# Sum over the whole time period ------------------------------------------

fls = list.files(inpath, pattern = '*.nc')

for (i in 1:length(fls)){
  
  f = nc_open(paste0(inpath, fls[i]))
  varsize <- f[["var"]][["net_infiltration"]][["size"]]
  
  df = array(0, dim = c(varsize[1], varsize[2]))
  for(j in 1:varsize[3]) df = df + ncvar_get(f, 'net_infiltration', start = c(1, 1, j), count = c(varsize[1], varsize[2], 1))
  df = df*0.0254 #meters
  
  size = round(ncvar_get(f, 'x')[2] - ncvar_get(f, 'x')[1])
  xll = ncvar_get(f, 'x')[1] - size/2
  yll = ncvar_get(f, 'y')[1] - size/2
  
  fname = paste0(outpath, '/', strsplit(fls[i], '_')[[1]][1], "_sum_tot.asc")
  save_ArcGRID(round(df, 6), fname, xll, yll, size, -9999)
  
  nc_close(f)
}
