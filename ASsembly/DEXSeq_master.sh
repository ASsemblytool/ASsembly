#!/bin/bash

# DEXseq master script
#
# Runs two pythonscripts to consecutively "flatten" the gtf-file and convert bam-files into count files.
# Then runs Rscript to perform analysis and parse results into a table of genes and p-values.

# Load defaults, most importantly regex names for control and treatment. Rest is mostly overwritten in the next lines, don't think that will cause problems.
source ASsembly_defaults.conf

DATA_PATH=$1 #Specify where to look for necessary data; BAM-files, gtf/gff-file
NCORES=$2 # Specify number of cores to use in R-part of scripti
GTF=$3 # Specify gtf-file
OUTPUT_PATH=$4/DEXSeq_output # Path where output should be put
mkdir -p $OUTPUT_PATH

PY_PATH=$(Rscript DEXSeq_findPythonScripts.R | grep -Po '(?<=\").*(?=\")') # Retrieve location of python files 

# Make flattened GFF-file out of GTF-file.
python ${PY_PATH}/dexseq_prepare_annotation.py -r no $GTF ${OUTPUT_PATH}/DEXSeq_FlattenedFeatureFile.gff

# Make count files
mkdir -p $OUTPUT_PATH/HTSeqCount_files
for file in ${DATA_PATH}/bam_files/*.bam; do
fileName=$(basename $file| grep -Po '.*(?=\.)')
python ${PY_PATH}/dexseq_count.py -f bam ${OUTPUT_PATH}/DEXSeq_FlattenedFeatureFile.gff ${file}  ${OUTPUT_PATH}/HTSeqCount_files/${fileName}_count.txt
done

Rscript DEXSeqScript.R $DATA_PATH $NCORES $TREATMENT $CONTROL
