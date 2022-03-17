#Script to compute the sum of the variable
# computed by SWB2 over the whole period considered
#Author: Paolo Colombo

#Works with the variables:
# - 'net_infiltration'
# - 'gross_precipitation'
#It can be tried with other variables, it should work, just hasn't been tested

# Setup -------------------------------------------------------------------

#ENG: Write the path to the main folder
#IT: Inserisci il percorso alla cartella principale
setwd('c:/E-OBS-SWB2')

#ENG: Write the path to the folder where you want the sums to be saved
#IT: Inserisci il percorso alla cartella dove vuoi che le somme vengano salvate
outpath = "./Export/ASCII/Sums" #don't write the last /

#ENG: Enter the path to the folder where you saved the SWB2 output.
# In that folder, put only the files of the desired variables
#IT: Inserisci il percorso alla cartella dove hai salvato gli output di SWB2
# In quella cartella, metti solo i file della variabile desiderata
inpath = "./Data/SWB2_output/gross_precipitation" #with the last /

#ENG: Write the variable name
#IT: Inserisci il nome della variabile
variable = 'gross_precipitation'

# Library and custom function needed --------------------------------------

library(ncdf4)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('custom_functions.R')

# Sum over the whole time period ------------------------------------------

fls = list.files(inpath, pattern = '*.nc')

for (i in seq_len(length(fls))){
  
  f = nc_open(paste0(inpath, '/', fls[i]))
  varsize <- f[["var"]][[variable]][["size"]]

  df = array(0, dim = c(varsize[1], varsize[2]))
  for(j in seq_len(varsize[3])) df = df + ncvar_get(f, variable, start = c(1, 1, j), count = c(varsize[1], varsize[2], 1))
  df = df*0.0254 #meters
  
  size = round(ncvar_get(f, 'x')[2] - ncvar_get(f, 'x')[1])
  xll = ncvar_get(f, 'x')[1] - size/2
  yll = ncvar_get(f, 'y')[length(ncvar_get(f, 'y'))] - size/2
  
  fname = paste0(outpath, '/', strsplit(fls[i], '_')[[1]][1], "_", variable, "_sum_tot.asc")
  save_ArcGRID(round(df, 6), fname, xll, yll, size, -9999)
  
  nc_close(f)
}
