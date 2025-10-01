# Generating spatial genetic layer to incorporate in species distribution modelling

This directory contains code and input data for generating a spatial genetic layer from kriged (interpolated) ancestry coefficients across the seascape.

We will use DArTseq data for <i>Aipysurus laevis</i>. The spatial genetic layer we will generate is the extent of the northwest shelf. Input files were downloaded from the GitHub repository: [sea-snake-dart](https://github.com/a-lud/sea-snake-dart). Environmental layers were downloaded from MARSPEC (via R). The northwest shelf shapefile was obtained from Vinay Udyawer.

<br>

### Prepare sample list
First, we will prepare a sample list which includes individuals of <i>A. laevis</i> from the northwest shelf.

```r
library(tidyverse)

laevis_popmap <- 
  read.csv("./genetic_layer/laevis/sample_list/ALA-popmap-coords.csv",
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

### Not run ###
# write.table(laevis_vcf_keep,
#             file = "./genetic_layer/laevis/sample_list/ALA-nw.txt",
#             sep = ",",
#             row.names = FALSE,
#             quote = FALSE)
```

The written output will be used to subset the VCF file.

<br>

### Subset the VCF file
Next, we will subset the VCF file by individual using `vcftools`. In bash, we run:
```bash
vcftools --vcf ./genetic_layer/laevis/vcf_file/ALA-stringent.highQ.filtered.vcf --keep ./genetic_layer/laevis/sample-list/ALA-nw.txt --recode --stdout > ./genetic_layer/laevis/vcf_file/ALA-stringent.highQ.filtered.nw.keep.vcf
```
