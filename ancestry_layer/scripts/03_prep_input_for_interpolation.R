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
  require(tess3r)
  require(raster)
  require(sdmpredictors)
  require(sp)
  require(dismo)
  require(terra)
  
  
  
  # download the layers -- already saved in marspec dir
  mspec_layers <- load_layers(mspec_annual, datadir = "../../marspec/")
  
  # Introduce NW shelf boundary
  nw_shelf <- sf::st_read("../../nw_shapefile/NWShelf.shp", quiet = TRUE) %>% sf::st_transform(4326)
  
  mspec_spatrast <- raster::crop(mspec_layers, nw_shelf)
  mspec_spatrast <- terra::rast(mspec_spatrast) # global; convert RasterBrick to SpatRaster
  mspec_shelf <- terra::mask(mspec_spatrast, mask = terra::vect(nw_shelf))
  #plot(mspec_shelf[[14]]) # 14 = sea surface temp (annual mean)
  
  mspec_sh_pcs <- RStoolbox::rasterPCA(mspec_shelf, spca = TRUE)
  mspec_sh_stack <- raster::stack(mspec_sh_pcs$map)
  #mspec_sh_stack
  
  mspec_sh_krig <- (mspec_sh_stack[[1]]) # environmental PC1, nw shape
  y_krig_sh_raster <- projectRaster(mspec_sh_krig, crs = "epsg:3112")
  


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
  colnames(species_keep_nw_samples) <- "sample"

  species_keep_nw_coords <- species_nw %>% 
    # arrange `laevis_keep_nw_coords` as in samples in VCF file
    arrange(id_clean, species_keep_nw_samples$sample) %>% 
    # select only long and lat cols
    select("longitude","latitude")
  colnames(species_keep_nw_coords) <- c("x","y") # rename cols
  
  # project
  species_keep_nw_proj <- 
    st_as_sf(species_keep_nw_coords, coords = c("x","y"), crs = "epsg:4326")
           
  species_keep_nw_proj <- 
    st_transform(species_keep_nw_proj, crs = 4326)
  
  message("Prep complete.")
  message("Objects generated: 1) imputed dosage, 2) projected coordinates, 3) extent for interpolation.")
  message("Proceed to TESS + Kriging - interpolation of ancestry coefficients across the seascape.")
  
  species_tess <- tess_ktest(gen = species_keep_nw_dosimp,
                             coords = species_keep_nw_proj,
                             Kvals = 1:length(unique(species_nw$pop)),
                             ploidy = 2,
                             K_selection = "auto")
  
  # Get TESS object and best K from results
  species_tessobj <- species_tess$tess3_obj
  species_bestK <- species_tess[["K"]]
  message(paste0("Automatic selection: Best K = ", species_bestK))
  
  # Get matrix of ancestry coefficients
  message("1...")
  species_qmat <- tess3r::qmatrix(tess3 = species_tessobj, K = species_bestK)
  
  # sample coordinates
  message("2...")
  x_species_keep_nw_coords <- st_as_sf(species_keep_nw_coords, coords = c("x", "y"), crs = 4326)
  x_proj_species_keep_nw_coords <- st_transform(x_species_keep_nw_coords, crs = 3112) # EPSG:3112; GDA94 Geoscience Australia projection 
  
  message("3...")
  species_krig_sh_admix <- 
    tess_krig(qmat = species_qmat, 
              coords = x_proj_species_keep_nw_coords, 
              grid = y_krig_sh_raster)
  
  message("4...")
  species_map_sh_admix <- 
    algatr::tess_ggplot(species_krig_sh_admix,
                        plot_method = "maxQ",
                        plot_axes = TRUE,
                        coords = x_proj_species_keep_nw_coords)
  
  plot(species_map_sh_admix)
  
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

laevis_tess <- 
  algatr::tess_ktest(gen = laevis_dos_imp, 
                     coords = laevis_nw_proj, 
                     Kvals = 1:7, 
                     ploidy = 2, 
                     K_selection = "auto")
# ^note that K changes every iteration ...

# Get TESS object and best K from results (i.e., laevis_tess)
laevis_tessobj <- laevis_tess$tess3_obj
laevis_bestK <- laevis_tess[["K"]]

# Get matrix of ancestry coefficients
laevis_qmat <- 
  tess3r::qmatrix(tess3 = laevis_tessobj, 
                  K = laevis_bestK)

# Prepare for kriging

# grid
mspec_sh_krig <- (mspec_sh_stack[[1]]) # environmental PC1, nw shape
y_krig_sh_raster <- raster::projectRaster(mspec_sh_krig, crs = "epsg:3112")

# sample coordinates
x <- sf::st_as_sf(laevis_nw_coords, coords = c("x", "y"), crs = 4326)
x_proj <- sf::st_transform(x, crs = 3112) # EPSG:3112; GDA94 Geoscience Australia projection 

z_krig_sh_admix <- 
  algatr::tess_krig(qmat = laevis_qmat, 
                    coords = x_proj, 
                    grid = y_krig_sh_raster)

z_map_sh_admix <- 
  algatr::tess_ggplot(z_krig_sh_admix,
                      plot_method = "maxQ",
                      plot_axes = TRUE,
                      coords = x_proj)
plot(z_map_sh_admix)
