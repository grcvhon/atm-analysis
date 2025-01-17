### MaxEnt modelling framework ###
### 2025-01-16 ###

### Install packages ---------------------------------------
# install.packages("sf")
# install.packages("leaflet")
# install.packages("readr")
# install.packages("janitor")
# install.packages("dplyr")
# install.packages("mapview")
# install.packages("spatstat")
# install.packages("tidyverse")
# install.packages("raster")
# install.packages("dismo")
# install.packages("lubridate")
# install.packages("SDMtune")
# install.packages("readxl")
# install.packages("stars")
# install.packages("lwgeom")

### Load packages -------------------------------------------
library(sf)
library(leaflet)
library(readr)
library(janitor)
library(dplyr)
library(mapview)
library(spatstat)
library(tidyverse)
library(raster)
library(dismo)
library(lubridate)
library(SDMtune)
library(readxl)
library(terra)
library(stars)
library(lwgeom)

### Input model extent ------------------------------------------

# Load Australia EEZ shapefile
eez <- st_read("misc/from-LM/australia-eez/Australia_EEZ.shp")
mapview(eez)

# Load environmental raster
bathymetry <- rast("../Data/Data from old models/Predictor variables/bathymetry.asc")
mapview(bathymetry)

# Clip the shapefile based on the raster
sf_clipped <- intersect(vect(eez), vect(ext(bathymetry))) # clip vectorised shapefile to the extent of the vectorised raster
mapview(sf_clipped) # shapefile clipped to the extent of the raster file (here used bathymetry.asc)


# Input occurrence data
dat <- read_excel("data/ATM_2023_0715-running-master.xlsx")
dat
