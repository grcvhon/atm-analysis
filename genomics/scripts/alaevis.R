# alaevis.R
# 06 June 2025 - combined from other script files

# Prepare sample list general - A laevis NW shelf samples ----
library(tidyverse)
# Load csv
laevis_popmap <- 
  read.csv("./genomics/ALA-popmap-coords.csv", 
           sep = ",", 
           header = TRUE)

laevis_nw <- laevis_popmap %>% 
  # Remove AL401 and AL404 - dropped samples
  subset(id!="AL401") %>% 
  subset(id!="AL404") %>% 
  # Keep only samples west of GoC
  subset(pop!="North_QLD") %>% 
  subset(pop!="Gulf_of_Carpentaria") %>% 
  subset(pop!="New_Caledonia")

laevis_nw <- laevis_nw %>%
  # unite first 3 column values separated with "-"
  unite("id_clean", species:targetid, remove = FALSE, sep = "-") %>% 
  # arrange columns as desired
  select("id_clean","pop","locality","longitude","latitude")
laevis_nw

# create keep list for vcf
laevis_vcf_keep <- laevis_nw %>% 
  select("id_clean")
# remove header from keep list
colnames(laevis_vcf_keep) <- NULL

# already written to file
#write.table(laevis_vcf_keep,
#            file = "./genomics/ALA-nw.txt",
#            sep = ",",
#            row.names = FALSE,
#            quote = FALSE)

# Genetic data processing - loading A laevis DArTseq VCF ----

# after filtering VCF via "Filter by individual" (i.e., use of keep list)
# load vcf file
laevis_vcf <- 
  vcfR::read.vcfR("./genomics/vcf_files/ALA-stringent.highQ.filtered.nw.keep.vcf", 
                  verbose = TRUE)

laevis_dosage <- 
  wingen::vcf_to_dosage(laevis_vcf)

laevis_dos_imp <- 
  algatr::simple_impute(laevis_dosage, FUN = median)

# Prepare sample list for algatr analyses ----

# obtain list of samples from VCF
laevis_samples <- as.data.frame(colnames(laevis_vcf@gt)) # get sample names from ala_vcf
laevis_samples <- as.data.frame(laevis_samples[-1,]) # remove "FORMAT" colname
colnames(laevis_samples) <- "sample" # rename "FORMAT" as "sample"
laevis_samples

laevis_nw_coords <- laevis_nw %>%
  # arrange `laevis_nw_coords` as in samples in vcf file
  arrange(id_clean, laevis_samples$sample) %>% 
  # select only long and lat cols
  select("longitude","latitude")
colnames(laevis_nw_coords) <- c("x","y") # rename cols

laevis_nw_proj <- sf::st_as_sf(laevis_nw_coords, coords = c("x", "y"), crs = "epsg:4326")
laevis_nw_proj <- sf::st_transform(laevis_nw_proj, crs = 4326)

# MARSPEC - get marine layer datasets ----

library(sdmpredictors)
library(raster)
library(sp)
library(dismo)

# MARSPEC data sets
# list names of MARSPEC layers (annual data) 
mspec_names <- list_layers("MARSPEC", monthly = FALSE)$name 
# all annual layers
mspec_annual <- list_layers("MARSPEC", monthly = FALSE)$layer_code 
# download the layers -- already saved in marspec dir
# mspec_layers <- load_layers(mspec_annual, datadir = "./genomics/marspec/") 

# Introduce NW shelf boundary ----
nw_shelf <- sf::st_read("./data/shapefiles/nw-shelf/NWShelf.shp", quiet = TRUE) %>% sf::st_transform(4326)
mapview::mapview(nw_shelf)

# Limit to NW shelf boundary ----
mspec_spatrast <- raster::crop(mspec_layers, nw_shelf)
mspec_spatrast <- terra::rast(mspec_spatrast) # global; convert RasterBrick to SpatRaster
mspec_shelf <- terra::mask(mspec_spatrast, mask = terra::vect(nw_shelf))
plot(mspec_shelf[[14]]) # 14 = sea surface temp (annual mean)
points(laevis_nw_coords, pch = 19)

mspec_sh_pcs <- RStoolbox::rasterPCA(mspec_shelf, spca = TRUE)
mspec_sh_stack <- raster::stack(mspec_sh_pcs$map)
mspec_sh_stack

# TESS + Kriging - interpolate ancestry coefficients across the land(sea)scape ----

# Notes:
# > In tess_ktest(), best K value changes every iteration (K_selection = "auto").
# > When K = 2, both Ks are represented on the plot. When K = 3, only K1 shows up. 

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

t <- raster::raster(z_krig_sh_admix$K1) # extent is different; does not represent coordinates
t <- terra::rast(t)
df_t <- as.data.frame(t, xy = TRUE)
write.csv(df_t, "./alaevis_tess_overlay.csv")

