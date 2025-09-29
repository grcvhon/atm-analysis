#######################
# generate keep files
#######################

# Files: 
# -- `*-popmap-coords.csv` are from github.com/a-lud/sea-snake-dart
# -- removed samples (not based on location) are from `samples-remove.txt` 
#    in github.com/a-lud/sea-snake-dart

# set work dir
setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/ancestry_layer/")

# load package
library(tidyverse)

##############
# laevis
##############

# Load csv
laevis_popmap <- 
  read.csv("./sample-lists/ALA-popmap-coords.csv", 
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
write.table(laevis_vcf_keep,
            file = "./sample-lists/ALA-nw.txt",
            sep = ",",
            row.names = FALSE,
            quote = FALSE)

##############
# major
##############

# Load csv
major_popmap <- 
  read.csv("./sample-lists/HMA-popmap-coords.csv",
           sep = ",", header = TRUE)
major_popmap
unique(major_popmap$pop) # determine non-nw pops

major_nw <- major_popmap %>% 
  # Remove samples based from "samples-remove.txt"
  subset(id!="VU46") %>% 
  subset(id!="SS170818_01") %>% 
  subset(id!="SS170814_02") %>% 
  subset(id!="SS170812_01") %>% 
  subset(id!="WAM-R154980") %>% 
  # remove non-nw samples
  subset(pop!="Gulf_of_Carpentaria") %>%  
  subset(pop!="Fraser_Island") %>%  
  subset(pop!="New_Caledonia")

major_nw

major_nw <- major_nw %>%
  # unite first 3 column values separated with "-"
  unite("id_clean", species:targetid, remove = FALSE, sep = "-") %>% 
  # arrange columns as desired
  # specify dplyr select 
  # otherwise conflicts with select function in raster package
  dplyr::select("id_clean","pop","locality","longitude","latitude")

major_nw

# create keep list for vcf
major_vcf_keep <- major_nw %>% 
  dplyr::select("id_clean")
# remove header from keep list
colnames(major_vcf_keep) <- NULL

# Not run - already written
write.table(major_vcf_keep,
            file = "./sample-lists/HMA-nw.txt",
            sep = ",",
            row.names = FALSE,
            quote = FALSE)

##############
# stokesii
##############

# Load csv
stokesii_popmap <- 
  read.csv("./sample-lists/HST-popmap-coords.csv",
           sep = ",", header = TRUE)
stokesii_popmap
unique(stokesii_popmap$pop) # determine non-nw pops

stokesii_nw <- stokesii_popmap %>% 
  # Remove samples based from "samples-remove.txt"
  subset(id!="KLS0360") %>% 
  subset(id!="KLS0634") %>% 
  subset(id!="KLS0886") %>% 
  subset(id!="KLS1679") %>% 
  subset(id!="KLS1677") %>% 
  subset(id!="KLS0660") %>% 
  subset(id!="KLS0941") %>% 
  subset(id!="KLS0897") %>% 
  subset(id!="KLS0896") %>% 
  subset(id!="KLS1204") %>% 
  subset(id!="As012") %>% 
  subset(id!="As010") %>% 
  subset(id!="KLS1688") %>% 
  # remove non-nw samples
  subset(pop!="Gulf_of_Carpentaria") %>%  
  subset(pop!="North_QLD")

stokesii_nw

stokesii_nw <- stokesii_nw %>%
  # unite first 3 column values separated with "-"
  unite("id_clean", species:targetid, remove = FALSE, sep = "-") %>% 
  # arrange columns as desired
  # specify dplyr select 
  # otherwise conflicts with select function in raster package
  dplyr::select("id_clean","pop","locality","longitude","latitude")

stokesii_nw

# create keep list for vcf
stokesii_vcf_keep <- stokesii_nw %>% 
  dplyr::select("id_clean")
# remove header from keep list
colnames(stokesii_vcf_keep) <- NULL

# Not run - already written
write.table(stokesii_vcf_keep,
            file = "./sample-lists/HST-nw.txt",
            sep = ",",
            row.names = FALSE,
            quote = FALSE)

# Next step:
# Use the keep list output to subset VCF files using vcfR