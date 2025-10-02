#!/usr/bin/env bash

#SBATCH --job-name=sbs_speed
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

#Rscript sbs_bearing_script.R
Rscript sbs_speed_script.R
