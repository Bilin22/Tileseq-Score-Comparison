#!/bin/bash

for i in {151..174}; do
    cat ${i}_*_L001_R1_001.fastq.gz ${i}_*_L002_R1_001.fastq.gz ${i}_*_L003_R1_001.fastq.gz ${i}_*_L004_R1_001.fastq.gz > /home/rothlab/bnong/tileseq/CBS/FASTQ/CBS_FASTQ_combined/818${i}_R1_joint.fastq.gz
    cat ${i}_*_L001_R2_001.fastq.gz ${i}_*_L002_R2_001.fastq.gz ${i}_*_L003_R2_001.fastq.gz ${i}_*_L004_R2_001.fastq.gz > /home/rothlab/bnong/tileseq/CBS/FASTQ/CBS_FASTQ_combined/818${i}_R2_joint.fastq.gz
done


for i in {181..192}; do
    cat ${i}_*_L001_R1_001.fastq.gz ${i}_*_L002_R1_001.fastq.gz ${i}_*_L003_R1_001.fastq.gz ${i}_*_L004_R1_001.fastq.gz > /home/rothlab/bnong/tileseq/CBS/FASTQ/CBS_FASTQ_combined/818${i}_R1_joint.fastq.gz
    cat ${i}_*_L001_R2_001.fastq.gz ${i}_*_L002_R2_001.fastq.gz ${i}_*_L003_R2_001.fastq.gz ${i}_*_L004_R2_001.fastq.gz > /home/rothlab/bnong/tileseq/CBS/FASTQ/CBS_FASTQ_combined/818${i}_R2_joint.fastq.gz
done