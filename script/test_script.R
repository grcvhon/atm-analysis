## MaxEnt modelling framework ##
## 2025-01-16 ##

#### Install packages ####
# install.packages("sf")
# install.packages("leaflet")
# install.packages("readr)
# install.packages("janitor")
# install.packages("dplyr")
# install.packages("mapview")
# install.packages("spatstat")
# install.packages("tidyverse")
# install.packages("raster")
# install.packages("dismo")
# install.packages("lubridate")
# install.packages("SDMtune")


#### Load packages ####
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

# Load shapefile

# Input occurrence data
occdata <- read_csv("data/ATM_2023_0715-running-master.csv")
occdata
