#!/bin/bash
cd /home/rothlab/bnong/tileseq/MTHFR/MTHFR_A222V_split_2023-07-13-15-58-18
touch check.txt
output=/home/rothlab/bnong/tileseq/MTHFR/MTHFR_A222V_split_2023-07-13-15-58-18/check.txt
cd sam_files
for SAMPLE in $(ls *R1_joint.sam|cut -d_ -f1); do
  LOG=$(ls -t ${SAMPLE}-3*.log)
  LAST=$(tail -n 1 ${LOG})
  echo ${LAST}>> ${output}
  R1LINES=$(wc -l ${SAMPLE}_R1_joint.sam|cut -f1 -d' ')
  R2LINES=$(wc -l ${SAMPLE}_R2_joint.sam|cut -f1 -d' ')
  if (( R1LINES == R2LINES )); then
    echo "Sample $SAMPLE OK" >> ${output}
  else
    echo "Sample $SAMPLE failed: $R1LINES lines in R1; $R2LINES lines in R2" >> ${output}
  fi
  echo "" >> ${output}
done
echo "DONE"