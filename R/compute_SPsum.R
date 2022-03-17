#Script to compute the sum of the variable in the stress periods provided
#Author: Paolo Colombo

#Works with the variables:
# - 'net_infiltration'
# - 'gross_precipitation'
#It can be tried with other variables, it should work, just hasn't been tested


# Setup -------------------------------------------------------------------

#Paths
#Direct path to the netCDF file you want to import
netCDFpath = 'C:/E-OBS-SWB2/Data/SWB2_output/0Impervious_net_infiltration.nc'
#Path where to save the ASCII files created
# Leave none if you don't want to create ASCII files
outpath = 'none'

#Specify the stress periods
SP1 = 90   #days, 01/01 - 30/03
SP2 = 76   #days, 01/04 - 12/06
SP3 = 92   #days, 13/06 - 15/09
SP4 = 107  #days, 16/09 - 31/12
SPs = c(SP1, SP2, SP3, SP4)
SPs = cumsum(SPs)

#Specify start and end year of the period
starty = 2014
endy = 2018

#Specify the variable
variable = 'net_infiltration'

# Library and custom function needed --------------------------------------

library(ncdf4)
library(tictoc)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('custom_functions.R')

# Sum over the specified stress periods -----------------------------------

SPsum_df = SPsum(path = netCDFpath, SPs = SPs, var = variable,
                 starty = starty, endy = endy,
                 outpath = outpath)
#It takes around 35 seconds for a 5 year, 660x338, 1day resolution netCDF

# Iterate over a list of files --------------------------------------------

#It doesn't work as it is
#Needs: changing outpath, automatic creation of new folders

inpath = ''
fls = list.files(inpath, pattern = '*.nc')
for (i in seq_len(length(fls))){
  SPsum(path = paste0(inpath, '/', fls[i]))
}

