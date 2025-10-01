# Generate spatial passage layer for species distribution modelling

This directory contains code and input data for generating a layer representing mean passage probability among spatially balanced points* across the northwest shelf. 

Mean passage probability was estimated based on the ocean current bearing^ which can produce asymmetrical routes between any two points. As such, pairwise mean passage probability was estimated and then visualised. The final output is a `.csv` file which can be used as input/predictor layer for species distribution modelling.

The code was written in R and executed using the University of Adelaide High Performance Computer (Phoenix). In this document, I present the code and visual output.

##
###### * Spatially balanced points plus manually selected points in Shark Bay and Exmouth Gulf to explicitly include such localities.<br>^ Mean passage probability estimates based on ocean current speed still needs to be executed.

<br>

### Download ocean current bearing and speed datasets
We will obtain our ocean current bearing and speed datasets from [BioOracle](https://www.bio-oracle.org/) via the R package `biooracler`.




### Run in Phoenix HPC
Bash script:
```bash
#!/usr/bin/env bash

#SBATCH --job-name=sbs
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH -o /hpcfs/users/a1235304/atm_passage/log/%x_%A.log

module load GDAL
module load GEOS
module load PROJ
module load UDUNITS
module load cURL
module load R

Rscript sbs_script.R
```
Run the bash script above using `sbatch <bash script>`.