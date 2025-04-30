#!/usr/bin/env bash
#SBATCH --job-name=06-pca
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --time=01:00:00
#SBATCH --mem=20GB
#SBATCH -o /hpcfs/users/a1235304/atm/slurm/06-pca/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vhon.garcia@adelaide.edu.au

PCA=/hpcfs/users/a1235304/atm/scripts/06-pca.py

source "/gpfs/apps/icl/software/Anaconda3/2024.06-1/etc/profile.d/conda.sh"
conda activate ipyrad

python ${PCA}