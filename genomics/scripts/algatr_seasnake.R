# algatr on sea snake data
# created: 09 May 2025
# by: Vhon Garcia

library(algatr)

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





