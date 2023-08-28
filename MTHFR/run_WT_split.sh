#!/bin/bash
submitjob.sh -n tileseqMut -c 8 -m 16G -t 4-00:00:00 -l MTHFR_WT_split.log -e MTHFR_WT_split.log -- \
tileseq_mut -p /home/rothlab/bnong/tileseq/MTHFR/parameters/MTHFR_WT_split.json \
-o /home/rothlab/bnong/tileseq/MTHFR/ -b galen1,galen2,galen3,galen4,galen70 -f /home/rothlab/bnong/tileseq/MTHFR/Data/MTHFR_FASTQ_joint \
--name MTHFR_WT_split --environment galen -c 6 --calibratePhredWT