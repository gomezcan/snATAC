#!/bin/bash

###################################
#######       Modules     #########
###################################

ml Bioinformatics  bcl2fastq2/2.20.0.422-oxq6lf3

###################################
###    Command lines to Run     ###
###################################

# variables/files input:
bcl_dic=$1     # Directory with bcl files: 

SampleIndex=$2 # SampleSheet.csv

 bcl2fastq --use-bases-mask=Y150n*,I8n*,Y16,Y150n* \
  --create-fastq-for-index-reads \
  --minimum-trimmed-read-length=8 \
  --mask-short-adapter-reads=8 \
  --ignore-missing-positions \
  --ignore-missing-controls \
  --ignore-missing-filter \
  --ignore-missing-bcls \
  -r 15 -w 15 -p 15\
  -R $bcl_dic \
  --output-dir=fastqs \
  --interop-dir=${bcl_dic}/InterOp \
  --sample-sheet=$SampleIndex
