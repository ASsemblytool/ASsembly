#!/bin/bash

# cufflinks script for AS recognition

# Load defaults, most importantly path to binaries. Rest is mostly overwritten in the next lines, don't think that will cause problems.
source ASsembly_defaults.conf

DATA_PATH=$1 # Path to data structure
BAM_PATH=$1/bam_files # Path to bam files
NCORES=$2 # Define number of cores to use
GTF=$3 # Define where to find the gtf-file
OUTPUT_PATH=$4 # Path where output should be put
tBams=$5 # comma-seperated string of treatment bam-files
cBams=$6 # comma-seperated string of control bam-files

mkdir -p ${OUTPUT_PATH}/Cufflinks_output/Cuffquant_output
mkdir -p ${OUTPUT_PATH}/Cufflinks_output/Cuffdiff_output
mkdir -p $OUTPUT_PATH/GSEA/Cufflinks

# Run cuffquant for every bam-file
for file in ${BAM_PATH}/*.bam; do
fileName=$(basename $file | grep -Po '.*(?=\.)')
${CUFFLINKS_PATH}/cuffquant -q -o ${OUTPUT_PATH}/Cufflinks_output/Cuffquant_output --library-type fr-firststrand -m $READLENGTH $GTF $file
mv ${OUTPUT_PATH}/Cufflinks_output/Cuffquant_output/abundances.cxb ${OUTPUT_PATH}/Cufflinks_output/Cuffquant_output/${fileName}.cxb
done

# add transcript start id (TSS_ID) to the gtf file
Rscript cuffdiff_gtf_attributes.R -i $GTF -o ${OUTPUT_PATH}/Cufflinks_output/gtf_wTSS.gtf

# Run cuffdiff
${CUFFLINKS_PATH}/cuffdiff -q -o ${OUTPUT_PATH}/Cufflinks_output/Cuffdiff_output -p $NCORES --library-type fr-firststrand -m $READLENGTH $OUTPUT_PATH/Cufflinks_output/gtf_wTSS.gtf $tBams $cBams

# Create overall toGSEA file
tail -n +2 ${OUTPUT_PATH}/Cufflinks_output/Cuffdiff_output/splicing.diff | awk '$7=="OK"' | cut -f2,12 > ${OUTPUT_PATH}/GSEA/Cufflinks/cufflinks_toGSEA.txt

