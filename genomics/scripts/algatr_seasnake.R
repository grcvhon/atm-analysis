# algatr on sea snake data
# created: 09 May 2025
# by: Vhon Garcia

library(algatr)
alazygatr_packages()

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
