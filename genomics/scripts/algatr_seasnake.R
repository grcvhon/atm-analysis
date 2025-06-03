# algatr on sea snake data
# created: 09 May 2025
# updated: 02 June 2025
# by: Vhon Garcia

# ---- preliminary ----
# load algatr library
library(algatr)
library(wingen)
library(raster)
library(terra)
library(ggplot2)
library(sf)

# load associated packages
alazygatr_packages()

# package example
load_algatr_example()
# loads `liz_vcf`,`liz_gendist`, `liz_coords`, `CA_env`

# wingen
# we reformat our dataframe of coordinates into sf coordinates
coords_longlat <- st_as_sf(liz_coords, coords = c("x","y"), crs = "+proj=longlat")

# the coordinates and raster can be projected to an equal area projection, in this case NAD83 / California Albers (EPSG 3310)
coords_proj <- st_transform(coords_longlat, crs = 3310)

# We'll also want the shape of California for subsequent plotting, 
# so we'll save one of the environmental PC layers as a SpatRaster object for this purpose and 
# reproject to the same coordinate reference system as the coordinates.
envlayer <- rast(CA_env[[1]])
# aggregate
envlayer <- aggregate(envlayer,5)
# reproject with same crs
envlayer <- project(envlayer, crs(coords_proj))

# Generate raster layer for sliding window
liz_lyr <- coords_to_raster(coords_proj, res = 50000, buffer = 5, plot = TRUE)

# get an idea of what the size of the cell and moving window look like using preview_gd()
sample_count <- preview_gd(liz_lyr, coords_proj, wdim = 3, fact = 0)

# visualise the sample count layer
ggplot_count(sample_count)

# run the moving window
wgd <- window_gd(liz_vcf,
                 coords_proj,
                 liz_lyr,
                 stat = "pi",
                 wdim = 3, fact = 0)

# visualise wingen results
ggplot_gd(wgd, bkg = envlayer) + ggtitle("Moving window pi")
ggplot_count(wgd) + ggtitle("Sample count")


kgd <- krig_gd(wgd, index = 1:2, liz_lyr, disagg_grd = 5)
ggplot_gd(kgd) + ggtitle("Kriged pi")


# read in AFO vcf
AFO_vcf <- read.vcfR("./genomics/vcf_files/AFO-reference.highQ.filtered.keep.vcf", verbose = TRUE)

# AFO_vcf
# ***** Object of Class vcfR *****
#   8 samples
# 155 CHROMs
# 6,249 variants
# Object size: 2 Mb
# 0 percent missing data
# *****        *****         *****

# obtain list of samples from VCF
AFO_keep <- as.data.frame(colnames(AFO_vcf@gt))
AFO_keep <- as.data.frame(AFO_keep[-1,])
colnames(AFO_keep) <- "sample"
AFO_keep

library(dplyr)

# load general sample sheet
sample_sheet <- read.delim("./genomics/sample-sheets/sample-sheet.csv", sep = ",", header = TRUE)
# filter general sample sheet to only contain samples in vcf
AFO_coords <- sample_sheet %>% filter(id_clean %in% AFO_keep$sample)
# reorder samples as in vcf
AFO_coords <- AFO_coords %>% arrange(id_clean, AFO_keep$sample)
# retain only longitude and latitude information
AFO_coords <- (AFO_coords[,8:9])
# move replicate coordinates by 0.0001; two snakes in one shot
AFO_coords[2,] <- AFO_coords[2,]+0.0001
# check that row 2 values have changed
AFO_coords
# rename colnames as x and y
colnames(AFO_coords) <- c("x","y")
# preview AFO_coords
AFO_coords


library(terra)
library(raster)
library(mapview)
library(sf)
library(maptools)
# read in extent (nw_shelf)
nw_shelf <- st_read("./data/shapefiles/nw-shelf/NWShelf.shp", quiet = TRUE) %>% st_transform(4326)
# read in environmental var (bathymetry)
bathymetry <- raster("./data/predictor-variables/bathymetry.asc")

# convert coordinates df to sf
# function
convert_2_sf <- function(df_name) {
  df_name %>% 
    st_as_sf(coords = c("x", "y"), crs = 4326) %>% 
    distinct(.keep_all = T)
}

# use function
AFO_points <- convert_2_sf(df_name = AFO_coords)
# visualise points on a map
mapview(AFO_points) + mapview(bathymetry, na.color=NA) # two points are the same x,y (two AFO in one shot = same GPS)


# TESS3
AFO_dosage <- vcf_to_dosage(AFO_vcf)
AFO_tess <- tess_ktest(AFO_dosage,
                       AFO_coords,
                       Kvals = 1:8,
                       ploidy = 2,
                       K_selection = "auto")
AFO_tess3obj <- AFO_tess$tess3_obj
AFO_K <- AFO_tess[["K"]]
AFO_qmat <- qmatrix(AFO_tess3obj, K = AFO_K)
AFO_qmat
tess_barplot(AFO_qmat)
tess_ggbarplot(AFO_qmat)

AFO_pcdists <- gen_dist(AFO_vcf, dist_type = "pc", npc_selection = "auto", criticalpoint = 0.9793)
gen_dist_hm(AFO_pcdists)


# MMRR
# Convert genetic data to matrix
Y <- as.matrix(AFO_pcdists)
# Extract values from our environmental raster
AFO_env <- raster::extract(bathymetry, AFO_points)
# Calculate environmental distances
X <- env_dist(AFO_env)
# Add geographic distance to X
X[["geodist"]] <- geo_dist(AFO_points)

# Run MMRR
AFO_resfull <- mmrr_run(Y, X, stdz = TRUE, nperm = 999, model = "full")
# Interpret results
AFO_resfull

# wingen - continuous mapping of genetic diversity using moving windows

AFO_coords # geographic coordinates of afo samples

# reformat dataframe of coordinates into sf coordinates
AFO_coords_longlat <- st_as_sf(AFO_coords, coords = c("x", "y"), crs = "+proj=longlat")
mapview(AFO_coords_longlat) # plot(AFO_coords_longlat)
AFO_coords_proj <- st_transform(AFO_coords_longlat, crs = 4326)
mapview(AFO_coords_proj) # plot(AFO_coords_proj)

# create raster from coordinates
AFO_lyr <- coords_to_raster(AFO_coords_proj, res = 1, buffer = 5, plot = TRUE)

# ---- Alaevis trial ---- 
# 2 June 2025

## Genetic data processing ----
library(algatr)
data_processing_packages()
library(purrr)
library(dplyr)
library(ggplot2)
library(here)

ala_vcf <- read.vcfR("./genomics/vcf_files/ALA-stringent.highQ.filtered.keep.vcf", verbose = TRUE)
ala_dosage <- vcf_to_dosage(ala_vcf)
ala_dos_imp <- simple_impute(ala_dosage, FUN = median)

## Environmental data processing ----
envirodata_packages()
library(terra)
library(raster)
library(RStoolbox)
library(ggplot2)
library(geodata)
library(viridis)
library(wingen)
library(tidyr)
library(tibble)

# obtain list of samples from VCF
ala_samples <- as.data.frame(colnames(ala_vcf@gt)) # get sample names from ala_vcf
ala_samples <- as.data.frame(ala_samples[-1,]) # remove "FORMAT" colname
colnames(ala_samples) <- "sample" # rename "FORMAT" as "sample"
ala_samples

ala_popmap <- read.delim("./genomics/ALA-popmap-coords.csv", header = TRUE, sep = ",") # take csv with long/lat info
ala_popmap <- ala_popmap %>% unite("id_clean", species:targetid, remove = FALSE, sep = "-") # unite first 3 column values separated with "-"
ala_popmap <- ala_popmap[-c(2,82),] # needed to single out AL401 and AL404
ala_coords <- ala_popmap %>% arrange(id_clean, ala_samples$sample) # arrange samples as in vcf
ala_coords <- ala_coords %>% select(longitude, latitude) # retain only long and lat columns
colnames(ala_coords) <- c("x","y") # rename cols

library(sf)
ala_coords <- st_as_sf(ala_coords, coords = c("x", "y"), crs = "+proj=longlat")
ala_proj <- st_transform(ala_coords, crs = 4326)

# download data from worldclim bound by ala_coords
ala_wclim <- get_worldclim(coords = ala_coords, res = 5, save_output = TRUE)
plot(ala_wclim[[1]])
points(ala_coords, pch = 19)

### MARSPEC ----
library(sdmpredictors)
marspec_annual_names <- list_layers("MARSPEC", monthly = FALSE)$name # as guide
marspec_annual <- list_layers("MARSPEC", monthly = FALSE)$layer_code # all annual layers
ala_marspec <- load_layers(marspec_annual)
#ala_marspec <- terra::rast(ala_marspec)
ala_marspec <- ala_marspec %>% raster::crop(ala_coords)
plot(ala_marspec[[14]])

## Detecting collinearity ----

# collinearity among environmental layers
# calculate the Pearson correlation coefficient for pairwise comparisons of environmental layers
cors_env <- check_env(ala_wclim) 
cors_env_m <- check_env(ala_marspec_crop)

# collinearity among extracted environmental variables
# determine collinearity using Pearson's correlation coefficients on extracted env vars at each sampling coordinate
# also generate a plot showing the pairwise correlations between environmental variables
check_result <- check_vals(ala_wclim, ala_coords) 
check_result_m <- check_vals(ala_marspec_crop, ala_coords)

# collinearity between distances
# determines collinearity between geographic and environmental distances
# does so by extracting values at sampling coordinates, 
# and then calculating geographic and environmental distances 
# and runs a Mantel test on resulting distances.
# Environmental distance = Euclidean distances
# Geographic distances = can be Euclidean, topographic, or resistance distance
check_results <- check_dists(ala_wclim, ala_coords)
check_results_m <- check_dists(ala_marspec_crop, ala_coords)
head(check_results$mantel_df)

## Raster PCA on environmental layers ----

# result of this function is a list containing the model information, 
# and a RasterBrick object containing multiple layers of PCA scores. 
# One can also use the nComp argument to only extract the top n PCs; 
# in many cases, the top three PCs may explain the majority of the variance of the data,
# and so only those will be considered for further analyses.

# using MARSPEC layers
env_pcs <- rasterPCA(ala_marspec, spca = TRUE)
envlayer <- rast(ala_marspec[[1]])
envlayer <- aggregate(envlayer, 5)
envlayer <- project(envlayer, crs(ala_proj))

# take a look at the results for the top 3 PCs
plots <- lapply(1:3, function(x) ggR(env_pcs$map, x, geom_raster = TRUE))
plots[[1]]
plots[[2]]
plots[[3]]

# We can also create a single composite raster plot with 3 PCs (each is assigned R, G, or B)
ggRGB(env_pcs$map, 1, 2, 3, stretch = "lin", q = 0, geom_raster = TRUE)


ala_lyr <- coords_to_raster(ala_proj, res = 0.1, buffer = 10, plot = TRUE)

ala_samp_count <- preview_gd(ala_lyr, ala_proj, wdim = 3, fact = 0)
ggplot_count(ala_samp_count)

ala_wgd <- window_gd(ala_vcf, ala_proj, ala_lyr, wdim = 3, fact = 0)

plot_gd(ala_wgd, bkg = envlayer)
plot_count(ala_wgd)


