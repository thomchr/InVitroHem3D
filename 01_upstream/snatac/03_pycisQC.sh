#!/bin/bash

#SBATCH --mail-type=END,FAIL 

source ~/.bashrc

conda activate scenicplus
dir=${EHT_DATA_DIR}/snATACseq
name=$1
pycistopic qc --fragments $dir/fragFiles/$name.tsv.gz --regions $dir/macs2_out/GRanges_clean.bed --tss $dir/cisTopic/qc/tss.bed --output $dir/cisTopic/qc/$name

