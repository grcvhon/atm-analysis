#######################################
# prepare input for TESS interpolation
#######################################

setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/ancestry_layer/vcf-files/nw/")

## Determine spatial extent of interpolation

# load marspec marine layer datasets
library(sdmpredictors)
library(raster)
library(sp)
library(dismo)

# MARSPEC data sets

# Not run
# list names of MARSPEC layers (annual data) 
# mspec_names <- list_layers("MARSPEC", monthly = FALSE)$name 

# all annual layers
# mspec_annual <- list_layers("MARSPEC", monthly = FALSE)$layer_code 

# download the layers -- already saved in marspec dir
mspec_layers <- load_layers(mspec_annual, datadir = "../../marspec/")

# Introduce NW shelf boundary
nw_shelf <- sf::st_read("../../nw_shapefile/NWShelf.shp", quiet = TRUE) %>% sf::st_transform(4326)
library(mapview)
mapview(nw_shelf)

mspec_spatrast <- raster::crop(mspec_layers, nw_shelf)
mspec_spatrast <- terra::rast(mspec_spatrast) # global; convert RasterBrick to SpatRaster
mspec_shelf <- terra::mask(mspec_spatrast, mask = terra::vect(nw_shelf))
plot(mspec_shelf[[14]]) # 14 = sea surface temp (annual mean)

mspec_sh_pcs <- RStoolbox::rasterPCA(mspec_shelf, spca = TRUE)
mspec_sh_stack <- raster::stack(mspec_sh_pcs$map)
mspec_sh_stack

## Prepare genetic species-specific input for interpolation

# Required objects from 01_generate_keep_files.R:
#   - laevis_nw
#   - major_nw
#   - stokesii_nw

prep_input <- function(keep_nw_vcf, species_nw){
  require(vcfR)
  require(wingen)
  require(algatr)
  require(dplyr)
  require(sf)

  laevis_nw
  major_nw
  stokesii_nw
  message("`laevis_nw`, `major_nw`, and `stokesii_nw` present in Environment.")

  species_keep_nw <- read.vcfR(keep_nw_vcf, verbose = FALSE)

  species_keep_nw_dosage <- vcf_to_dosage(species_keep_nw)

  species_keep_nw_dosimp <- simple_impute(species_keep_nw_dosage, FUN = median)

  # get samples names
  species_keep_nw_samples <- as.data.frame(colnames(species_keep_nw@gt))

  # remove "FORMAT" colname
  species_keep_nw_samples <- as.data.frame(species_keep_nw_samples[-1,])

  # rename as "sample"
  colnames(samples_keep_nw_samples) <- "sample"

  species_keep_nw_coords <- species_nw %>% 
    # arrange `laevis_keep_nw_coords` as in samples in VCF file
    arrange(id_clean, species_keep_nw_samples$sample) %>% 
  # select only long and lat cols
    select("longitude","latitude")
  colnames(species_keep_nw_coords) <- c("x","y") # rename cols
  
  # project
  species_keep_nw_proj <- st_as_sf(laevis_keep_nw_coords, coords = c("x","y"), crs = "epsg:4326")
           
  laevis_keep_nw_proj <- st_transform(laevis_keep_nw_proj, crs = 4326)

}

library(vcfR)
laevis_keep_nw <- 
  read.vcfR("./laevis_keep_nw.vcf", verbose = TRUE)

library(wingen)
laevis_keep_nw_dosage <- 
  vcf_to_dosage(laevis_keep_nw)

library(algatr)
laevis_keep_nw_dosimp <-
  simple_impute(laevis_keep_nw_dosage,
                FUN = median)

# prepare for algatr analyses

# obtain list of samples from VCF
laevis_keep_nw_samples <- as.data.frame(colnames(laevis_keep_nw@gt)) # get sample names
laevis_keep_nw_samples <- as.data.frame(laevis_keep_nw_samples[-1,]) # remove "FORMAT" colname
colnames(laevis_keep_nw_samples) <- "sample" # rename as "sample"
laevis_keep_nw_samples

library(dplyr)
laevis_keep_nw_coords <- laevis_nw %>% 
  # arrange `laevis_keep_nw_coords` as in samples in VCF file
  arrange(id_clean, laevis_keep_nw_samples$sample) %>% 
  # select only long and lat cols
  select("longitude","latitude")
colnames(laevis_keep_nw_coords) <- c("x","y") # rename cols

# project
library(sf)

laevis_keep_nw_proj <- 
  st_as_sf(laevis_keep_nw_coords, 
           coords = c("x","y"),
           crs = "epsg:4326")

laevis_keep_nw_proj <- 
  st_transform(laevis_keep_nw_proj, crs = 4326)