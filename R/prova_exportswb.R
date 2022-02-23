#Manage SWB2 netCDF output

# Setup -------------------------------------------------------------------

library(ncdf4)

#Project folder
setwd("C:/Users/paolo/Desktop/progetto E-OBS")

# Load outputs ------------------------------------------------------------

#Output folder
swbout = "./swb2_MODELMI/output/"
#List the netCDF files inside the folder
ls = list.files(path = swbout, pattern = "*.nc")



