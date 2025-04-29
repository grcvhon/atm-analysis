#!/usr/bin/env bash
#SBATCH --job-name=06-vcf2hdf5
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --time=1:00:00
#SBATCH --mem=30GB
#SBATCH -o /hpcfs/users/a1235304/atm/slurm/06-vcf2hdf5/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vhon.garcia@adelaide.edu.au

source "/gpfs/apps/icl/software/Anaconda3/2024.06-1/etc/profile.d/conda.sh"
conda activate ipyrad

python3

import ipyrad.analysis as ipa
