#!/usr/bin/env bash
#SBATCH --job-name=05-vcf-filter-highQ
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --time=1:00:00
#SBATCH --mem=60GB
#SBATCH -o /hpcfs/users/a1235304/atm/slurm/05-bcftools/%x_%j.log
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vhon.garcia@adelaide.edu.au

# Variables
DIR='/hpcfs/users/a1235304/atm'
DENOVO=$(find "${DIR}/results/ipyrad" -name '*-denovo.vcf.gz')
REFERENCE=$(find "${DIR}/results/ipyrad" -name '*-reference.vcf.gz')

source "/gpfs/apps/icl/software/Anaconda3/2024.06-1/etc/profile.d/conda.sh"
conda activate bcftools

for D in $DENOVO; do

    BN=$(basename "${D}" -denovo.vcf.gz)

    echo "${BN}"

    OUT="${DIR}/results/ipyrad/${BN}-denovo_outfiles"
    POP="${DIR}/data/popmaps/${BN}-popmap.tsv"

    bcftools view "${D}" |
        bcftools +setGT -- -t q -n . -i 'FMT/DP<5' |
        bcftools +fill-tags -- -t F_MISSING -S "${POP}" |
        bcftools filter -i 'F_MISSING<=0.1'| 
        bgzip > ${OUT}/${BN}-denovo.highQ.filtered.vcf.gz
    
    tabix ${OUT}/${BN}-denovo.highQ.filtered.vcf.gz

done

for R in $REFERENCE; do

    BNR=$(basename "${R}" -reference.vcf.gz)

    echo "${BNR}"

    OUTR="${DIR}/results/ipyrad/${BNR}-reference_outfiles"
    POPR="${DIR}/data/popmaps/${BNR}-popmap.tsv"

    bcftools view "${R}" |
        bcftools +setGT -- -t q -n . -i 'FMT/DP<5' |
        bcftools +fill-tags -- -t F_MISSING -S "${POPR}" |
        bcftools filter -i 'F_MISSING<=0.1'| 
        bgzip > ${OUTR}/${BNR}-reference.highQ.filtered.vcf.gz
    
    tabix ${OUTR}/${BNR}-reference.highQ.filtered.vcf.gz

done