library(sf)
library(tidyverse)

laevis_anc <- read.csv("./genomics/algatr_TESS/alaevis_tess_K2.csv", header = T)
laevis_anc <- laevis_anc[,-1]
laevis_anc_sf <- st_as_sf(laevis_anc, coords = c("x", "y"), crs = 4326)

laevis_flow <- read.csv("./genomics/scripts_passage/regular_seed1_50pts_17h53m10s/regular_seed1_50pts_17h53m10s.csv", header = T)
laevis_flow <- laevis_flow[,-1]
colnames(laevis_flow) <- c("x","y","p_pass")
laevis_flow_sf <- st_as_sf(laevis_flow, coords = c("x", "y"), crs = 4326)




