setwd("C:/Users/paolo/OneDrive/Documents/Dati")
library(raster)
r <- raster("fg_ens_mean_0.1deg_reg_2011-2021_v24.0e.nc")
install.packages("ncdf4")
library(ncdf4)
r <- nc_open("fg_ens_mean_0.1deg_reg_2011-2021_v24.0e.nc")
p <- nc_open("rr_ens_mean_0.1deg_reg_2011-2021_v24.0e.nc")


lat = ncvar_get(r, 'latitude')
lon = ncvar_get(ncf, 'longitude')

View(lat)
