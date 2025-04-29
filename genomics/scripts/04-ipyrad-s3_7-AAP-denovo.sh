#!/usr/bin/env bash
#SBATCH --job-name=AAP-denovo
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --time=1:00:00
#SBATCH --mem=60GB
#SBATCH -o /hpcfs/users/a1235304/atm/slurm/04-ipyrad-s37/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vhon.garcia@adelaide.edu.au

# Variables
DIR='/hpcfs/users/a1235304/atm'
cd "${DIR}/results/ipyrad" || exit 1

#source "/home/a1645424/hpcfs/micromamba/etc/profile.d/micromamba.sh"
#micromamba activate ipyrad

source "/gpfs/apps/icl/software/Anaconda3/2024.06-1/etc/profile.d/conda.sh"

conda activate ipyrad

ipyrad \
    -s 34567 \
    -p 'params-AAP-denovo.txt' \
    -c "${SLURM_CPUS_PER_TASK}"

conda deactivate

#micromamba deactivate
