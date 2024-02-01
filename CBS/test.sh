#!/bin/bash

for i in {151..174}; do
    cat ${i}_*_L001_R1_001.fastq.gz ${i}_*_L002_R1_001.fastq.gz ${i}_*_L003_R1_001.fastq.gz ${i}_*_L004_R1_001.fastq.gz > /home/rothlab/bnong/tileseq/CBS/FASTQ/CBS_FASTQ_combined/818${i}_R1_joint.fastq.gz
    cat ${i}_*_L001_R2_001.fastq.gz ${i}_*_L002_R2_001.fastq.gz ${i}_*_L003_R2_001.fastq.gz ${i}_*_L004_R2_001.fastq.gz > /home/rothlab/bnong/tileseq/CBS/FASTQ/CBS_FASTQ_combined/818${i}_R2_joint.fastq.gz
done


for i in {181..192}; do
    cat ${i}_*_L001_R1_001.fastq.gz ${i}_*_L002_R1_001.fastq.gz ${i}_*_L003_R1_001.fastq.gz ${i}_*_L004_R1_001.fastq.gz > /home/rothlab/bnong/tileseq/CBS/FASTQ/CBS_FASTQ_combined/818${i}_R1_joint.fastq.gz
    cat ${i}_*_L001_R2_001.fastq.gz ${i}_*_L002_R2_001.fastq.gz ${i}_*_L003_R2_001.fastq.gz ${i}_*_L004_R2_001.fastq.gz > /home/rothlab/bnong/tileseq/CBS/FASTQ/CBS_FASTQ_combined/818${i}_R2_joint.fastq.gz
done


# rerun the alignment files

#!/bin/bash

for i in {607004..607006}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

for i in {607008..607009}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

sabtch home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/607011_aln.sh

for i in {607030..607031}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

for i in {607033..607035}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

for i in {607016..607018}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

for i in {818037..818038}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

for i in {818048..818048}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

for i in {818020..818020}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

for i in {818128..818128}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

for i in {818060..818060}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done

for i in {818060..818023}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done


# outside for loop need to be re write

for i in {818060..818023}; do
    sbatch /home/rothlab/bnong/tileseq/CBS/CBS_2024-01-25-13-31-31/GALEN_jobs/${i}_aln.sh
done



