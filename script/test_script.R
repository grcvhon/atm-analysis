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

### Load package list ---------------------------------------
packages <- c("sf","leaflet","readr","janitor","dplyr",
              "mapview","spatstat","tidyverse","raster",
              "dismo","lubridate","SDMtune","readxl",
              "terra","stars","lwgeom")
invisible(lapply(packages, library, character.only = TRUE))

### Input model extent ------------------------------------------

# Load Northwest shelf bounding area
nw_shelf <- st_read("data/shapefiles/nw-shelf/NWShelf.shp")
mapview(nw_shelf)

# Load initial environmental predictors
env_init <- stack("data/predictor-variables/sal_mean.asc",
                  "data/predictor-variables/sal_amp.asc",
                  "data/predictor-variables/bathymetry.asc",
                  "data/predictor-variables/sst_mean.asc",
                  "data/predictor-variables/sst_amp.asc",
                  "data/predictor-variables/chlor_mean.asc",
                  "data/predictor-variables/DistToLand.asc",
                  "data/predictor-variables/DistToReef.asc",
                  "data/predictor-variables/DistToFW.asc")
env_init