#!/usr/bin/env bash

DIR='/hpcfs/users/a1235304/atm'
DENOVO=$(find "${DIR}/results/ipyrad" -name '*-denovo.vcf.gz')
REFERENCE=$(find "${DIR}/results/ipyrad" -name '*-reference.vcf.gz')

for D in $DENOVO; do
    BN=$(basename "${D}" -denovo.vcf.gz)

    echo "${BN}"

    OUT="${DIR}/results/ipyrad/${BN}-denovo_outfiles"
    POP="${DIR}/data/popmaps/${BN}-popmap.tsv"

    echo "${OUT}"
    echo "${POP}"
done