
install.packages("ncdf4")
library(ncdf4)

setwd("C:/Users/paolo/Desktop/progetto E-OBS/Dati/E-OBS/")

ls = list.files(pattern = '*.nc')

ncf = nc_open(ls[7]) #7 perch� file rr � il settimo nella lista
lat = ncvar_get(ncf, 'latitude')
lon = ncvar_get(ncf, 'longitude')

rr_mean = ncvar_get(ncf, 'rr', start = c(1, 1, 1), count = c(705, 465, 1))

row.names(t_mean) = lon
colnames(t_mean) = lat
