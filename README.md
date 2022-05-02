# E-OBS - SWB2
**Goal**: Management of E-OBS data to import them inside SWB2. Handling of SWB2 outputs.\
**Author**: Paolo Colombo - Politecnico di Milano\
**Contact**: paolocolombo1996@gmail.com

## What the codes do
- Handle E-OBS data
  - Creation of a class to handle E-OBS data: `EOBSobject`
- Support the creation of input files for SWB2 starting from E-OBS data:
  - Create netCDF files compatible with SWB2
  - Create ArcASCII GRID files compatible with SWB1 and SWB2
- Handle SWB2 netCDF output
  - Creation of a class to handle SWB2 output: `SWB2output`
- Calculate the total recharge by adding urban and irrigation components:
  - Creation of a class to perform the needed operations: `RechargeCalc`

## Repository organization
The repository is still undergoing massive reorganization to be functional and easily understandable.
It will be explained once it is done

## Useful links
SWB2 manual: https://pubs.er.usgs.gov/publication/tm6A59 \
E-OBS user guide: https://surfobs.climate.copernicus.eu/userguidance/use_ensembles.php \
E-OBS download page: https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php
