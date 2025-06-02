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
ala_coords

# download data from worldclim bound by ala_coords
ala_wclim <- get_worldclim(coords = ala_coords, res = 5, save_output = TRUE)
plot(ala_wclim)
points(ala_coords, pch = 19)
check_env(ala_wclim)

# ...collinearity issues next
