#!/usr/bin/env bash

# Note: this bash script is used to rename the processed fastq files (= undergone QC i.e., 02_qc.sh)
# however, the `02_qc.sh` step was done separately in a different Phoenix session/account (with a1645424)
# because the `kraken2` produces output with incorrect format with a1235304 Phoenix session/account
#
# so as a workaround, a1645424 ran the `02_qc.sh` script on fastq files of interest
# then provided processed files (i.e., QC'd) to a1235304
#
# a1235304 renamed these processed (QC'd) files to the desired filename format to proceed with the pipeline
#
# ***Troubleshooting ongoing as to why the same script, same `kraken2` version produces two different outputs***

PROCESSEDIR="/hpcfs/users/a1235304/atm/temp-vhon"
RENAMEDIR="/hpcfs/users/a1235304/atm/data/rename_processed_fastq"

mkdir -p "${RENAMEDIR}"

cd "${RENAMEDIR}" || exit 1

# Read 'sample-sheet.csv' line by line
while IFS=',' read -r -a LINE; do
  if [[ ! -f "${PROCESSEDIR}/${LINE[1]}.fastq.gz" ]]; then
    echo "File for ${LINE[1]} does not exist"
    exit 1
  fi
  echo -e "FROM: ${PROCESSEDIR}/${LINE[1]}.fastq.gz\tTO: ${LINE[2]}.fastq.gz"
  rsync -aP "${PROCESSEDIR}/${LINE[1]}.fastq.gz" .
  mv ./"${LINE[1]}.fastq.gz" ./"${LINE[2]}.fastq.gz"
done < <(tail -n +2 '/hpcfs/users/a1235304/atm/data/sample-sheet.csv')