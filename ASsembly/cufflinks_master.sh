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

# Run cuffquant for every bam-file
mkdir -p ${OUTPUT_PATH}/cuffquant_output

for file in ${BAM_PATH}/*.bam; do
fileName=$(basename $file | grep -Po '.*(?=\.)')
${CUFFLINKS_PATH}/cuffquant -q -o ${OUTPUT_PATH}/cuffquant_output --library-type fr-firststrand $GTF $file
mv ${OUTPUT_PATH}/cuffquant_output/abundances.cxb ${OUTPUT_PATH}/cuffquant_output/${fileName}.cxb
done

# Create output folder if necessary
mkdir -p ${OUTPUT_PATH}/cuffdiff_output

# add transcript start id (TSS_ID) to the gtf file
Rscript cuffdiff_gtf_attributes.R -i $GTF -o ${OUTPUT_PATH}/cuffdiff_output/gtf_wTSS.gtf

# Run cuffdiff
${CUFFLINKS_PATH}/cuffdiff -q -o ${OUTPUT_PATH}/cuffdiff_output -p $NCORES --library-type fr-firststrand $OUTPUT_PATH/cuffdiff_output/gtf_wTSS.gtf $tBams $cBams

# Create toGSEA file
tail -n +2 ${OUTPUT_PATH}/cuffdiff_output/splicing.diff | awk '$14=="yes"' | cut -f2,12 > ${OUTPUT_PATH}/toGSEA/cufflinks_toGSEA.txt
