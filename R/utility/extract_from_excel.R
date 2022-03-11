#Creazione di nuovi file utilizzabili per il calcolo della ricarica totale
#
# - Estrazione dei fogli e delle tabelle dal file "Ricarica_TesiPiccioli.xlsx"
# - Esportazione in CSV

# Setup -------------------------------------------------------------------

library("readxl")
setwd('C:/E-OBS-SWB2')

# Copy the sheet out of the big excel file --------------------------------
#And save it as single CSV files, more manageable

#Copia della prima tabella nel foglio "Calcolo aree" di Ricariche2
ind <- read.table('clipboard', sep = '\t', header = TRUE)
write.csv(ind, "./Data/Calcolo_ricarica_totale/indicatori.csv")

#Copia del foglio "Ricarica urbana" di Ricarica_TesiPiccioli
rurb <- read.table('clipboard', sep = '\t', header = TRUE)
write.csv(rurb, "./Data/Calcolo_ricarica_totale/ricarica_urbana.csv")

#Copia del foglio "Ricarica irrigua" di Ricarica_TesiPiccioli
#Il foglio ha due tabelle. Non comprendo quale delle due va 

#Copia foglio "Lam" di Ricarica_TesiPiccioli
lam <- read_excel('./Data/Calcolo_ricarica_totale/Ricarica_TesiPiccioli.xlsx', sheet = 7)
write.csv(lam, "./Data/Calcolo_ricarica_totale/lam.csv")



