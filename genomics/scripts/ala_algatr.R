# using algatr on alaevis DArTseq data
# created: 4 June 2025
# by: Vhon Garcia

library(tidyverse)

# Load csv
laevis_popmap <- read.csv("./genomics/ALA-popmap-coords.csv", 
                          sep = ",", header = TRUE)

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
write.table(laevis_vcf_keep,
            file = "./genomics/ALA-nw.txt",
            sep = ",",
            row.names = FALSE,
            quote = FALSE)

## Genetic data processing ----
library(vcfR)
laevis_vcf <- 
  read.vcfR("./genomics/vcf_files/ALA-stringent.highQ.filtered.nw.keep.vcf", 
            verbose = TRUE)

library(wingen)
laevis_dosage <- vcf_to_dosage(laevis_vcf)

library(algatr)
laevis_dos_imp <- simple_impute(laevis_dosage, FUN = median)

## Environmental data processing ----

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

library(sf)
laevis_nw_proj <- st_as_sf(laevis_nw_coords, coords = c("x", "y"), crs = "epsg:4326")
laevis_nw_proj <- st_transform(laevis_nw_proj, crs = 4326)

## --- laevis_nw_proj extext ---      ##
## xmin = 113.9903, xmax = 124.0431,  ##
## ymin = -26.39904, ymax = -12.1     ##

# get marine layers
library(sdmpredictors)
library(raster)
library(sp)
library(dismo)

# MARSPEC data sets
# list names of MARSPEC layers (annual data) 
mspec_names <- list_layers("MARSPEC", monthly = FALSE)$name 
# all annual layers
mspec_annual <- list_layers("MARSPEC", monthly = FALSE)$layer_code 
# download the layers
mspec_layers <- load_layers(mspec_annual) 
# plot raster layers
bound <- extent(c(xmin = 111.9903, 
                  xmax = 126.0431,
                  ymin = -28.39904, 
                  ymax = -11.1))
mspec_raster <- crop(mspec_layers, bound)
plot(mspec_raster)
points(laevis_nw_proj)

library(RStoolbox)
mspec_pcs <- rasterPCA(mspec_raster, spca = TRUE)
mspec_pcs
# this mspec_stack, same class as `CA_env` in accompanying example
mspec_stack <- raster::stack(mspec_pcs$map) 

# TESS
laevis_tess <- tess_ktest(laevis_dosage, laevis_nw_proj, Kvals = 1:7, ploidy = 2, K_selection = "auto")
# ^note that K changes every iteration ... 

# Get TESS object and best K from results (...laevis_tess)
laevis_tessobj <- laevis_tess$tess3_obj
laevis_bestK <- laevis_tess[["K"]]

# Get qmat of ancestry coeff
laevis_qmat <- tess3r::qmatrix(laevis_tessobj, K = laevis_bestK)

# tess3 krig raster
mspec_krig <- (mspec_stack[[1]]) # no aggregate factor

laevis_krig_admix <- tess_krig(laevis_qmat, laevis_nw_proj, mspec_krig)
