#!/bin/bash

# DEXseq master script
#
# Runs two pythonscripts to consecutively "flatten" the gtf-file and convert bam-files into count files.
# Then runs Rscript to perform analysis and parse results into a table of genes and p-values.

# Load defaults, most importantly regex names for control and treatment. Rest is mostly overwritten in the next lines, don't think that will cause problems.
source ASsembly_defaults.conf

DATA_PATH=$1 #Specify where to look for necessary data; BAM-files, gtf/gff-file
NCORES=$2 # Specify number of cores to use in R-part of script
GTF=$3 # Specify gtf-file
TREATMENT=$4
CONTROL=$5
INT=$6 # Path where output should be put
mkdir -p $INT/DEXSeq_output

PY_PATH=$(Rscript DEXSeq_findPythonScripts.R | grep -Po '(?<=\").*(?=\")') # Retrieve location of python files 

# Make flattened GFF-file out of GTF-file.
GFF="${INT}/DEXSeq_output/DEXSeq_FlattenedFeatureFile.gff"
python ${PY_PATH}/dexseq_prepare_annotation.py -r no $GTF $GFF

# Make count files
mkdir -p $INT/DEXSeq_output/HTSeqCount_files
for file in $DATA_PATH/bam_files/*.bam; do
fileName=$(basename $file| grep -Po '.*(?=\.)')
python ${PY_PATH}/dexseq_count.py -f bam $INT/DEXSeq_output/DEXSeq_FlattenedFeatureFile.gff ${file}  $INT/DEXSeq_output/HTSeqCount_files/${fileName}_count.txt
done

Rscript DEXSeqScript.R $INT $GFF $NCORES $TREATMENT $CONTROL
