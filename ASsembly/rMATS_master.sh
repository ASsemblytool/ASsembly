#!/bin/bash

# rMATS script
#


# Load defaults, most importantly path to binaries. Rest is mostly overwritten in the next lines, don't think that will cause problems.
source ASsembly_defaults.conf

DATA_PATH=$1 # Specify where to look for necessary data; BAM-files, GTF-file
GTF=$3 # Path to gtf-file
OUTPUT_PATH=$4 # Path where output should be put
mkdir -p $OUTPUT_PATH/rmats_output
tBams=$5 # comma-seperated string of treatment bam-files
cBams=$6 # comma-seperated string of control bam-files
READLENGTH=$7 #read length
PAIRED=$8
if [ $PAIRED=="TRUE" ]; then
PAIRED=paired
else
PAIRED=single
fi

BAM_PATH=$1/bam_files # Path to bam files.

# Run MATS
python ${RMATS_PATH}/RNASeq-MATS.py -b1 $tBams -b2 $cBams -gtf $GTF -o $OUTPUT_PATH/rmats_output -t $PAIRED -len $READLENGTH -analysis U -novelSS 1

# Extract AS events from files
for file in ${OUTPUT_PATH}/rmats_output/MATS_output/*.txt ; do
asType=$(basename $file | grep -o -P '^[^\.]*(?=\.)') # Extract type of AS from file name
asc=$(wc -l $file) # Count number of AS cases for this category
if [ "$asc" > "1" ]; then
tail -n +2 $file | cut -f3,19 | awk -v var="$asType" '{print $0"	"var}' | sort -g -k2,2 > ${OUTPUT_PATH}/toGSEA/${asType}_toGSEA.txt # add cases to file
fi
done
