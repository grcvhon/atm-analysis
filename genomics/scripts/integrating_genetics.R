# A Landscape Genomic Analysis Toolkit in R (algatr)
# based on Chambers, E.A., Bishop, A.P., & Wang, I.J. (2023). Individual-based landscape genomics for conservation: 
# An analysis pipeline. Molecular Ecology Resources.https://doi.org/10.1111/1755-0998.13884.

# Install package
# devtools::install_github("TheWangLab/algatr")

# Load package
library(algatr)

# Install all packages for algatr
# This installs all the algatr dependencies
# to perform analyses
alazygatr_packages()

load_algatr_example()

plot(CA_env, col = turbo(100), axes = FALSE)
env <- scaleRGB(CA_env)
plotRGB(env, r = 1, g = 2, b = 3)
points(liz_coords, pch = 19)

# For your own dataset, you need:
# > Sampling coordinates (longitude, latitude [x,y] order)
# > Genetic data in vcf file format
# > Environmental data layers
# !! Individuals in coordinates data file should be !!
# !! ordered the same way as in the genetic data file !!

# Above, we have loaded environmental data and individual coordinates data.
# The genetic data file gets loaded using read.vcfR()

### TESS3 in algatr ####

# Requirements: genotype dosage matrix, coordinates for samples, environmental layers
# use `vcf_to_dosage()` function to convert a vcf to a dosage matrix

# Convert vcf to genotype dosage matrix
liz_dosage <- vcf_to_dosage(liz_vcf)

krig_raster <- raster::aggregate(CA_env[[1]], fact = 6)
Qraster <- terra::plot(CA_env[[1]], col = mako(100), axes = FALSE)
terra::plot(krig_raster, col = mako(100), axes = FALSE)
points(liz_coords)

tess3_result <- tess_ktest(liz_dosage, # genotype dosage matrix
                           liz_coords, # coordinates object
                           Kvals = 1:10, # evaluate between K values 1 to 10
                           ploidy = 2,
                           K_selection = "auto")

tess3_obj <- tess3_result$tess3_obj
bestK <- tess3_result[["K"]]

qmat <- qmatrix(tess3_obj, K = bestK)
qmat
tess_barplot(qmat)
tess_ggbarplot(qmat)

coords_proj <- sf::st_as_sf(liz_coords, coords = c("x", "y"), crs = 4326)
coords_proj <- sf::st_transform(coords_proj, crs = 3310)
krig_raster <- raster::projectRaster(krig_raster, crs = "epsg:3310")
krgi_admix <- tess_krig(qmat, # qmatrix with ancestry coefficient values for each individual and each K value
                        coords_proj, # projected coordinates into coordinate system 
                        krig_raster) # kriging raster reprojected in the same coordinate system
krig_admix <- tess_krig(qmat, 
                        coords_proj, 
                        krig_raster)
terra::plot(krig_admix)

### Genetic distances ####
gen_dist_packages()
library(cowplot)

# Calculate genetic distances
pc_dists <- gen_dist(liz_vcf, dist_type = "pc", npc_selection = "auto", criticalpoint = 2.0234)
gen_dist_hm(pc_dists)
pc_dists

devtools::install("popmaps-1.03", dependency = FALSE)
library(popmaps)

