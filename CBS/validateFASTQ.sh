#!/bin/bash

printReadIDS() {
  FASTQ=$1
  zcat "$FASTQ"|awk '{if (NR % 4 == 1) {print $1}}'
}

SAMPLES=$(for FQ in *fastq.gz; do echo "${FQ%%_*}"; done|sort|uniq)
for SAMPLE in $SAMPLES; do
  printf "FASTQ for sample %s : " $SAMPLE
  FASTQ1=$(ls ${SAMPLE}_*_R1*.fastq.gz)
  FASTQ2=$(ls ${SAMPLE}_*_R2*.fastq.gz)
  ISDIFF=$(diff -q <(printReadIDS "$FASTQ1") <(printReadIDS "$FASTQ2"))
  if [[ "$ISDIFF" == Files*differ ]]; then
    echo "❌"
  else
    echo "✅"
  fi
done


#!/bin/bash

printReadIDS() {
  FASTQ=$1
  zcat "$FASTQ"|awk '{if (NR % 4 == 1) {print $1}}'
}

SAMPLES=$(for FQ in *fastq.gz; do echo "${FQ%%_*}"; done|sort|uniq)
for SAMPLE in $SAMPLES; do
  printf "FASTQ for sample %s : " $SAMPLE
  FASTQ1=$(ls ${SAMPLE}_R1*.fastq.gz)
  FASTQ2=$(ls ${SAMPLE}_R2*.fastq.gz)
  ISDIFF=$(diff -q <(printReadIDS "$FASTQ1") <(printReadIDS "$FASTQ2"))
  if [[ "$ISDIFF" == Files*differ ]]; then
    echo "❌"
  else
    echo "✅"
  fi
done
