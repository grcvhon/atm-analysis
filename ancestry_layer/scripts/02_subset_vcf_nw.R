##################################
# subset updated* VCF using vcfR
##################################

subset_vcf <- function(vcf, keep_list, save_as_vcf, save_as_gzvcf){
  
  require(vcfR)
  
  # set work dir
  setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/ancestry_layer/")
  
  # read vcf (supply updated vcf)
  readvcf <- read.vcfR(vcf, verbose = FALSE)
  message("FOR COMPARISON: VCF **BEFORE** subsetting:")
  print(readvcf)
  
  # read keep list
  readkeep <- readLines(keep_list)
  message(paste0("`", keep_list, "` has ", length(readkeep)," samples."))
  
  # subset
  readvcf@gt <- readvcf@gt[, c("FORMAT", readkeep)]
  
  # set write dir  
  setwd("C:/Users/a1235304/Dropbox/Short-nosed and Leaf-scaled sea snake TSSC/atm-analysis/ancestry_layer/vcf-files/nw/")
  
  # write output
  write.vcf(readvcf, save_as_vcf)
  write.vcf(readvcf, save_as_gzvcf)
  message("~ ~ ~ Subsetting summary ~ ~ ~")
  print(readvcf)
  message("Output written to atm-analysis/ancestry_layer/vcf-files/nw/")
  
}

# laevis
subset_vcf(vcf = "./vcf-files/updated/ALA-stringent.vcf.gz",
           keep_list = "./sample-lists/ALA-nw.txt",
           save_as_vcf = "./laevis_keep_nw.vcf",
           save_as_gzvcf = "./laevis_keep_nw.vcf.gz")

# major
subset_vcf(vcf = "./vcf-files/updated/HMA-stringent.vcf.gz",
           keep_list = "./sample-lists/HMA-nw.txt",
           save_as_vcf = "./major_keep_nw.vcf",
           save_as_gzvcf = "./major_keep_nw.vcf.gz")

# stokesii
subset_vcf(vcf = "./vcf-files/updated/HST-stringent.vcf.gz",
           keep_list = "./sample-lists/HST-nw.txt",
           save_as_vcf = "./stokesii_keep_nw.vcf",
           save_as_gzvcf = "./stokesii_keep_nw.vcf.gz")

