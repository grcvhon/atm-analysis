#!/usr/bin/env bash

#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=72:00:00
#SBATCH --mem=32GB
#SBATCH -o /hpcfs/users/a1235304/atm/slurm/%x_%j.log

# Notification configuration
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=vhon.garcia@adelaide.edu.au

RESSTORE='/uofaresstor/sanders_lab/sequencing-datasets/radseq'
FQDIR="/hpcfs/users/${USER}/atm/data/fastq"

mkdir -p "${FQDIR}"

cd "${FQDIR}" || exit 1

# Read 'sample-sheet.csv' line by line
while IFS=',' read -r -a LINE; do
  if [[ ! -f "${RESSTORE}/DaRT-${LINE[0]}/${LINE[1]}.FASTQ.gz" ]]; then
    echo "File for ${LINE[2]} does not exist in ${LINE[0]}"
    exit 1
  fi
  echo -e "FROM: ${RESSTORE}/DaRT-${LINE[0]}/${LINE[1]}.FASTQ.gz\tTO: ${LINE[2]}.FASTQ.gz" 
  rsync -aP "${RESSTORE}/DaRT-${LINE[0]}/${LINE[1]}.FASTQ.gz" .
  rename -v "${LINE[1]}" "${LINE[2]}" "${LINE[1]}.FASTQ.gz"
done < <(tail -n +2 '/hpcfs/users/a1235304/atm/data/sample-sheet.csv')