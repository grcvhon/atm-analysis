# tutorial: gdistance package for randomise least-cost paths
# https://www.alexbaecher.com/post/connectivity-script/

library(gdistance)
library(tidyverse)
library(rgeos)
library(elevatr)
library(ggplot2)
library(tigris)
library(spocc)
library(raster)
library(viridis)
library(ggthemes)

# tigris package for USA census data - salamander
states <- tigris::states()
se <- states %>% subset(REGION == "3")
se
TN_NC <- se %>% subset(NAME %in% c("Tennessee", "North Carolina")) # subset TN and NC data

Pj <- occ(query = "Plethodon jordani",                 # JRCS scientific name
          from = "gbif",                               # limiting query to *the first* 1000 records
          limit=1000,                                  # limiting query to *the first* 1000 records
          has_coords = T)                              # limiting those 1000 records to those that have geo-referenced data

Pj_sp <- Pj$gbif$data$Plethodon_jordani %>%            # Grabbing the Darwin-core data from the spocc object
  dplyr::select(longitude,                             # Keep locations and year, discard the rest
                latitude,
                year) %>%                   
  dplyr::filter(year > 2000) %>%                       # Filter records to only those after year 2000
  filter(!duplicated(round(longitude, 2),              # Remove duplicate records using rounded decimals (this removes points very near to one-another)       
                     round(latitude, 2)) == TRUE) %>%  # >> See notes below about ^^
  dplyr::mutate(lon = scale(longitude),                # Remove points far outside the cluster of occurrences
                lat = scale(latitude)) %>%             # >> See notes below about ^^
  dplyr::filter(!abs(lon)>2) %>%
  dplyr::filter(!abs(lat)>2) %>%
  dplyr::select(longitude,
                latitude) %>%
  SpatialPoints(proj4string = crs(se)) 

# thinking about it... the resistance values may not need optimisation
# because, given site A and B and their bearings from each other
# there would be a calculated value based on the locations of the sites and their relative pairwise bearings
# taking from Wilcox et al 2023
# Resistance was then calculated as the absolute difference in the current deviation from the deviation of the 
# median bearing to the desired bearings and standardized by dividing by the greatest potential deviance (e.g., resistance=|deviation—90|/90).
# also
# Values for resistance were calculated based on the deviation of current direction from the bearing angle to the destination site. 
# A deviation value of 0° represents no resistance (current direction=bearing angle) and a deviation value of 180° represents maximum resistance, 
# which would be standardized for resistance through division of the maximum potential deviation (e.g., 180°).
# 
# 
# Vinay Udyawer - That would make it easier to calculate
# 
# and take out the guesswork as well
# then we can standardise the values across the extent of analyis
# 
# 
# Vinay Udyawer - Although since we want to predict across the full model space.. the resistance layer needs to be generalized across the full model space
# 
# yes, the values of the bearing in the cells will help with that
# for example, even if A is 45 degrees to B, the cells in between may have different bearing values
# so res.AtoB is not necessarily = 0
# 
# Vinay Udyawer - Right but how do we define A and B? Will it just be occurrence locations?
# 
# A and B will be the occurrence locations, yes but it will be represented by sample locations (= populations)
# 
# Vinay Udyawer - That sounds reasonable
# 
# say Pilbara (that has x samples), then Exmouth (y number of samples) and so on
# so taking bearing between Pilbara and Exmouth, and Exmouth to Pilbara
# if Pilbara is -135 degrees to Exmouth, the path of least resistance will be the path with the least deviation to the bearing. 
# such path is determined by the bearing values in the raster grid. so it can result to a map showing 
# multiple corridors - some better (least/lesser resistance), some not as much (some resistance)