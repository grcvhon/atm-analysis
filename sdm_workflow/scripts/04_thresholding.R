# sdm - 04 thresholding

# determine suitable and non-suitable areas from the predicted distribution 
# based on a threshold value

## setwd: sdm directory
#setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/sdm")

# required packages
library(terra)
library(sf)
library(tidyverse)


# short-nosed sea snake

threshold_snss <- results_list$short_nosed_opt_mod.rds$thresholds[3,2] # cloglog of max training sens and spec

# coastal thresholding
short_nosed_PA <- short_nosed_predict_list$coastal$prediction_mean # coastal mean prediction
short_nosed_PA[short_nosed_PA < threshold_snss] <- NA # below threshold
short_nosed_PA[short_nosed_PA >= threshold_snss] <- 1 # above threshold
plot(short_nosed_PA)

short_nosed_PA_poly <- as.polygons(short_nosed_PA)
plot(short_nosed_PA_poly) # SpatVector
short_nosed_PA_area <- st_area(st_as_sf(short_nosed_PA_poly)) # area in m^2
short_nosed_PA_area <- short_nosed_PA_area/1e+06 # convert to km^2
# SNSS (coastal) = 30,527.75 square kilometres


# offshore thresholding 
short_nosed_offPA <- short_nosed_predict_list$offshore$prediction_mean # offshore mean prediction
short_nosed_offPA[short_nosed_offPA < threshold_snss] <- NA # below threshold
short_nosed_offPA[short_nosed_offPA >= threshold_snss] <- 1 # above threshold
plot(short_nosed_offPA)

short_nosed_offPA_poly <- as.polygons(short_nosed_offPA)
plot(short_nosed_offPA_poly) # SpatVector
short_nosed_offPA_area <- st_area(st_as_sf(short_nosed_offPA_poly)) # area in m^2
short_nosed_offPA_area <- short_nosed_offPA_area/1e+06 # convert to km^2
# SNSS (offshore) = 4,009.76 square kilometres

ggplot() + 
  geom_spatraster(data = terra::rast(nw_shelf)) +
  geom_spatraster(data = short_nosed_PA, na.color = NA) +
  geom_spatraster(data = short_nosed_offPA, na.color = NA) +
  #labs(title = paste0("Mean spatial prediction for Short-nosed sea snake (", last_part, ")")) +
  theme_bw()



# leaf-scaled sea snake

threshold_lsss <- results_list$leaf_scaled_opt_mod.rds$thresholds[3,2]

# coastal
leaf_scaled_coastPA <- leaf_scaled_predict_list$coastal$prediction_mean
leaf_scaled_coastPA[leaf_scaled_coastPA < threshold_lsss] <- NA # below threshold
leaf_scaled_coastPA[leaf_scaled_coastPA >= threshold_lsss] <- 1 # above threshold
plot(leaf_scaled_coastPA)

leaf_scaled_coastPA_poly <- as.polygons(leaf_scaled_coastPA)
plot(leaf_scaled_coastPA_poly) #SpatVector
leaf_scaled_coastPA_area <- st_area(st_as_sf(leaf_scaled_coastPA_poly)) # area in m^2
leaf_scaled_coastPA_area <- leaf_scaled_coastPA_area/1e+06 # convert to km^2
# LSSS (coastal) = 8,549.944 square kilometres


# offshore
leaf_scaled_offPA <- leaf_scaled_predict_list$offshore$prediction_mean
leaf_scaled_offPA[leaf_scaled_offPA < threshold_lsss] <- NA # below threshold
leaf_scaled_offPA[leaf_scaled_offPA >= threshold_lsss] <- 1 # above threshold
plot(leaf_scaled_offPA)


leaf_scaled_offPA_poly <- as.polygons(leaf_scaled_offPA)
plot(leaf_scaled_offPA_poly) # SpatVector
leaf_scaled_offPA_area <- st_area(st_as_sf(leaf_scaled_offPA_poly)) # area in m^2
leaf_scaled_offPA_area <- leaf_scaled_offPA_area/1e+06 # convert to km^2
# LSSS (offshore) = 1,393.315 square kilometres





#######################################################################################
### MPA overlap
#######################################################################################

# Collaborative Australian Protected Areas Database (CAPAD) 2024 - Marine
# Source: https://fed.dcceew.gov.au/datasets/0b6e7b6c48a64a3a82c225fa48aee13d_1/explore?location=-40.711703%2C115.622969%2C3

# read in capad2024 marine
d <- sf::st_read("../data/capad2024_marine/capad2024_marine.shp", quiet = TRUE) %>% 
  st_transform(4326)
# convert to SpatVector
d_vect <- terra::vect(d)

# read in nw shelf region
nw_shelf <- st_read("../data/shapefiles/nw-shelf/NWShelf.shp", quiet = TRUE) %>% 
  st_transform(4326)
# convert to SpatVector
nw_vect <- terra::vect(nw_shelf)

# retain MPA regions within model extent
d_crop <- terra::crop(d_vect, nw_vect)

# convert to sf
d_sf <- st_as_sf(d_crop)

# view
mapview(d_sf %>% group_by(EPBC) %>% summarise(geometry)) + 
  mapview(nw_shelf, alpha.regions = 0)


# calculate EPBC-specific areas

# cth
d_cth <- d_sf %>% filter(EPBC == "Commonwealth")
d_cth <- st_union(d_cth) # unite cth features
d_cth_area <- st_area(d_cth) # area of cth parks
d_cth_area <- d_cth_area/1e+06

d_cth_sf <- st_as_sf(d_cth)


# indigenous
d_ind <- d_sf %>% filter(EPBC == "Indigenous")
d_ind <- st_union(d_ind)
d_ind_area <- st_area(d_ind)
d_ind_area <- d_ind_area/1e+06

d_ind_sf <- st_as_sf(d_ind)


# state  
d_ste <- d_sf %>% filter(EPBC == "State")
d_ste_valid <- st_make_valid(d_ste)
d_ste_buffer <- st_buffer(d_ste_valid, 0)
d_ste_union <- st_union(d_ste_buffer)
d_ste_u_val <- st_make_valid(d_ste_union)
d_ste_u_area <- st_area(d_ste_u_val)/1e+6

d_ste_sf <- st_as_sf(d_ste_u_val)


# prepare species-specific prediction sf
lsss_off_sf <- st_as_sf(leaf_scaled_offPA_poly)
lsss_cst_sf <- st_as_sf(leaf_scaled_coastPA_poly)

snss_off_sf <- st_as_sf(short_nosed_offPA_poly) 
snss_cst_sf <- st_as_sf(short_nosed_PA_poly)

# intersections
leaf_scaled_int_off_cth <- st_intersection(lsss_off_sf, d_cth_sf)
leaf_scaled_int_cst_cth <- st_intersection(lsss_cst_sf, d_cth_sf)

leaf_scaled_int_off_ste <- st_intersection(lsss_off_sf, d_ste_sf)
leaf_scaled_int_cst_ste <- st_intersection(lsss_cst_sf, d_ste_sf)

leaf_scaled_int_off_ind <- st_intersection(lsss_off_sf, d_ind_sf)
leaf_scaled_int_cst_ind <- st_intersection(lsss_cst_sf, d_ind_sf)

short_nosed_int_off_cth <- st_intersection(snss_off_sf, d_cth_sf)
short_nosed_int_cst_cth <- st_intersection(snss_cst_sf, d_cth_sf)

short_nosed_int_off_ste <- st_intersection(snss_off_sf, d_ste_sf)
short_nosed_int_cst_ste <- st_intersection(snss_cst_sf, d_ste_sf)

short_nosed_int_off_ind <- st_intersection(snss_off_sf, d_ind_sf)
short_nosed_int_cst_ind <- st_intersection(snss_cst_sf, d_ind_sf)

# preview intersection
map(intersection = short_nosed_int_cst_cth,
    epbc = d_cth_sf, 
    prediction = snss_cst_sf)

# function: preview intersection
map <- function(intersection,epbc,prediction){
  mapview(intersection, col.regions = "chartreuse") + 
  mapview(epbc, alpha.regions = 0, color = "yellow", lwd = 2) + 
  mapview(prediction, col.regions = "darkolivegreen")
  }


intersections <- 
  list(
    LS_offshore_commonwealth = leaf_scaled_int_off_cth,
    LS_coastal_commonwealth = leaf_scaled_int_cst_cth,
    LS_offshore_state = leaf_scaled_int_off_ste,
    LS_coastal_state = leaf_scaled_int_cst_ste,
    LS_offshore_indigenous = leaf_scaled_int_off_ind,
    LS_coastal_indigenous = leaf_scaled_int_cst_ind,
    SS_offshore_commonwealth = short_nosed_int_off_cth,
    SS_coastal_commonwealth = short_nosed_int_cst_cth,
    SS_offshore_state = short_nosed_int_off_ste,
    SS_coastal_state = short_nosed_int_cst_ste,
    SS_offshore_indigenous = short_nosed_int_off_ind,
    SS_coastal_indigenous = short_nosed_int_cst_ind
    )

intersection_area <- lapply(intersections, FUN = function(x) st_area(x) / 1e6)
intersection_area


ggplot() + 
  geom_sf(data = nw_shelf, fill = "cadetblue1") + 
  geom_sf(data = lsss_off_sf, fill = "pink", col = NA) +
  geom_sf(data = lsss_cst_sf, fill = "pink", col = NA) +
  geom_sf(data = d_ste_sf, col = "cyan3", fill = NA) +
  geom_sf(data = d_cth_sf, col = "cyan4", fill = NA) +
  geom_sf(data = leaf_scaled_int_cst_ste, fill = "red", col = NA) +
  geom_sf(data = leaf_scaled_int_off_cth, fill = "red", col = NA) +
  annotation_scale(mapping = aes(location = "br")) +
  #coord_sf(xlim = c(110, 116), ylim = c(-27.5,-21), expand = FALSE) +
  theme_bw()

ggplot() + 
  geom_sf(data = nw_shelf, fill = "cadetblue1") + 
  geom_sf(data = snss_off_sf, fill = "pink", col = NA) +
  geom_sf(data = snss_cst_sf, fill = "pink", col = NA) +
  geom_sf(data = d_ste_sf, col = "cyan3", fill = NA) +
  geom_sf(data = d_cth_sf, col = "cyan4", fill = NA) +
  geom_sf(data = short_nosed_int_cst_ste, fill = "red", col = NA) +
  geom_sf(data = short_nosed_int_cst_cth, fill = "red", col = NA) +
  geom_sf(data = short_nosed_int_off_cth, fill = "red", col = NA) +
  annotation_scale(mapping = aes(location = "br")) +
  #coord_sf(xlim = c(110, 116), ylim = c(-27.5,-21), expand = FALSE) +
  theme_bw()

#######################################################################################
### Fisheries overlap
#######################################################################################
# Shark Bay Prawn Limited Entry Fishery (DPIRD-062) (CMP)
# Source: https://services.slip.wa.gov.au/public/rest/services/SLIP_Public_Services/DPIRD_Fisheries_Guide/MapServer/77

sb_fish <- sf::st_read("../data/fisheries/Shark Bay Prawn Limited Entry Fishery (DPIRD-062) (CMP).shp", quiet = TRUE) %>% 
  st_transform(4326)
sb_fish_vect <- terra::vect(sb_fish)
sb_fish_sf <- st_as_sf(sb_fish_vect)

# The fishery
sb_fish_1 <- sb_fish_sf[1,]
# closure areas
sb_fish_2 <- sb_fish_sf[2,] 
sb_fish_3 <- sb_fish_sf[3,]
sb_fish_4 <- sb_fish_sf[4,]

# remove closure areas from The fishery
sb_fish_only <- terra::erase(terra::vect(sb_fish_1), terra::vect(sb_fish_2))
sb_fish_only <- terra::erase(sb_fish_only, terra::vect(sb_fish_3))
sb_fish_only <- terra::erase(sb_fish_only, terra::vect(sb_fish_4))
sb_fish_only <- st_as_sf(sb_fish_only)


# Exmouth Gulf Prawn Managed Fishery (DPIRD-062) (CMP)
# Used here from Vinay
# Also available at: https://services.slip.wa.gov.au/public/rest/services/SLIP_Public_Services/DPIRD_Fisheries_Guide/MapServer/61
eg_fish <- sf::st_read("../data/fisheries/Exmouth_trawl_ground.geojson") %>% 
  st_transform(4326)
eg_fish_vect <- terra::vect(eg_fish)
eg_fish_sf <- st_as_sf(eg_fish_vect)

eg_fish_1 <- eg_fish_sf[1,] # The fishery
eg_fish_2 <- eg_fish_sf[2,] # nursery
eg_fish_3 <- eg_fish_sf[3,] # gear trial
eg_fish_4 <- eg_fish_sf[4,] # port area

eg_fish_only <- terra::erase(terra::vect(eg_fish_1), terra::vect(eg_fish_2))
eg_fish_only <- terra::erase(eg_fish_only, terra::vect(eg_fish_3))
eg_fish_only <- terra::erase(eg_fish_only, terra::vect(eg_fish_4))
eg_fish_only <- st_as_sf(eg_fish_only)

sb_fish_only_area <- st_area(sb_fish_only)/1e+6
eg_fish_only_area <- st_area(eg_fish_only)/1e+6

# intersection
leaf_scaled_int_cst_sbfish <- st_intersection(lsss_cst_sf, sb_fish_only)
leaf_scaled_int_cst_sbfish_area <- st_area(leaf_scaled_int_cst_sbfish)/1e+6

leaf_scaled_int_cst_egfish <- st_intersection(lsss_cst_sf, eg_fish_only)
leaf_scaled_int_cst_egfish_area <- st_area(leaf_scaled_int_cst_egfish)/1e+6

short_nosed_int_cst_sbfish <- st_intersection(snss_cst_sf, sb_fish_only)
short_nosed_int_cst_sbfish_area <- st_area(short_nosed_int_cst_sbfish)/1e+6

short_nosed_int_cst_egfish <- st_intersection(snss_cst_sf, eg_fish_only)
short_nosed_int_cst_egfish_area <- st_area(short_nosed_int_cst_egfish)/1e+6

# preview intersection
map(intersection = short_nosed_int_cst_egfish,
    epbc = eg_fish_only, 
    prediction = snss_cst_sf)

library(rworldmap)
library(rworldxtra)
world <- getMap(resolution = "high")
world_sf <- st_as_sf(world)

wa <- sf::st_read("../data/STE_2021_AUST_SHP_GDA94/STE_2021_AUST_GDA94.shp") %>% 
  st_transform(4326)
wa <- wa[5,]
wa <- st_as_sf(wa)


# snss sb
ggplot() + 
  geom_sf(data = wa, fill = "lightyellow") + 
  geom_sf(data = snss_cst_sf, fill = "pink", col = NA) +
  geom_sf(data = sb_fish_only, col = "darkblue", fill = NA) +
  geom_sf(data = short_nosed_int_cst_sbfish, fill = "red", col = NA) +
  annotation_scale(mapping = aes(location = "bl")) +
  coord_sf(xlim = c(112, 114.5), ylim = c(-26.9,-23.3), expand = FALSE) +
  theme_pubclean() +
  theme(panel.background = element_rect(fill = "lightcyan1"),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank()) 

# lsss sb
ggplot() + 
  geom_sf(data = wa, fill = "lightyellow") + 
  geom_sf(data = lsss_cst_sf, fill = "pink", col = NA) +
  geom_sf(data = sb_fish_only, col = "darkblue", fill = NA) +
  geom_sf(data = leaf_scaled_int_cst_sbfish, fill = "red", col = NA) +
  annotation_scale(mapping = aes(location = "bl")) +
  coord_sf(xlim = c(112, 114.5), ylim = c(-26.9,-23.3), expand = FALSE) +
  theme_pubclean() +
  theme(panel.background = element_rect(fill = "lightcyan1"),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank()) 

# snss eg
ggplot() + 
  geom_sf(data = wa, fill = "lightyellow") + 
  geom_sf(data = snss_cst_sf, fill = "pink", col = NA) +
  geom_sf(data = eg_fish_only, col = "darkblue", fill = NA) +
  geom_sf(data = short_nosed_int_cst_egfish, fill = "red", col = NA) +
  annotation_scale(mapping = aes(location = "br")) +
  coord_sf(xlim = c(114, 114.8), ylim = c(-22.6,-21.5), expand = FALSE) +
  theme_pubclean() +
  theme(panel.background = element_rect(fill = "lightcyan1"),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) 

# lsss eg
ggplot() + 
  geom_sf(data = wa, fill = "lightyellow") + 
  geom_sf(data = lsss_cst_sf, fill = "pink", col = NA) +
  geom_sf(data = eg_fish_only, col = "darkblue", fill = NA) +
  geom_sf(data = leaf_scaled_int_cst_egfish, fill = "red", col = NA) +
  annotation_scale(mapping = aes(location = "br")) +
  coord_sf(xlim = c(114, 114.8), ylim = c(-22.6,-21.5), expand = FALSE) +
  theme_pubclean() +
  theme(panel.background = element_rect(fill = "lightcyan1"),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))