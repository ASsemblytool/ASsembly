#!/bin/bash

# ASpli master script
#
# Runs Rscript to perform analysis and parse results into a table of genes and p-values,
# together with the event that happens (ES, IR, 3SS or 5SS).

# Load defaults, most impoartanly regex names for control and treatment.
source ASsembly_defaults.conf

DATA_PATH=$1 #Specify where to look for necessary data; BAM-files, gtf/gff-file
GTF=$2 #Specify where the GTF file can be found
INT=$3 #Path where output should be stored
REGEX_CONTROL=$4 #Specify the regex of the control files
REGEX_TREATMENT=$5 #Specify the regex of the treatment files
mkdir -p $INT/ASpli_output

Rscript ASpliScript.R $DATA_PATH $GTF $REGEX_TREATMENT $REGEX_CONTROL $READLENGTH $MAXINTRONLENGTH $INT/ASpli_output $NCORES
