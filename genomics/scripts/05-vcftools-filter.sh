#!/usr/bin/env bash
#SBATCH --job-name=05-vcftools-filter
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --time=1:00:00
#SBATCH --mem=80GB
#SBATCH -o /hpcfs/users/a1235304/atm/slurm/05-vcftools/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vhon.garcia@adelaide.edu.au

# load conda
source "/gpfs/apps/icl/software/Anaconda3/2024.06-1/etc/profile.d/conda.sh"
conda activate vcftools

# variables
DIR=/hpcfs/users/a1235304/atm/results/
IPYRAD=ipyrad/AFO-reference_outfiles/
VCF_IN=${DIR}${IPYRAD}AFO-reference.vcf.gz

# filters
MAC=3
MISS=0.9
QUAL=30
MIN_DEPTH=5
MAX_DEPTH=30

VCF_OUT=${DIR}${VCFTOOLS}AFO-reference.q${QUAL}.min${MIN_DEPTH}.max${MAX_DEPTH}.miss${MISS}.mac${MAC}.vcf.gz

vcftools --gzvcf ${VCF_IN} --remove-indels --mac ${MAC} --max-missing ${MISS} --minQ ${QUAL} --min-meanDP ${MIN_DEPTH} --max-meanDP ${MAX_DEPTH} --minDP ${MIN_DEPTH} --maxDP ${MAX_DEPTH} --recode --stdout | gzip -c > ${VCF_OUT}