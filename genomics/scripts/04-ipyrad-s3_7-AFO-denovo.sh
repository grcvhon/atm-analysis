#!/usr/bin/env bash
#SBATCH --job-name=AFO-denovo
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --time=4:00:00
#SBATCH --mem=60GB
#SBATCH -o /hpcfs/users/a1235304/atm/slurm/ipyrad-s37/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vhon.garcia@adelaide.edu.au

# Variables
DIR='/hpcfs/users/a1235304/atm'
cd "${DIR}/results/ipyrad" || exit 1

#source "/home/a1645424/hpcfs/micromamba/etc/profile.d/micromamba.sh"
#micromamba activate ipyrad
conda activate ipyrad
ipyrad \
    -s 34567 \
    -p 'params-AFO-denovo.txt' \
    -c "${SLURM_CPUS_PER_TASK}"
conda deactivate
#micromamba deactivate
