#!/bin/bash

# ASsembly script
#
# Assures input is valid, applies any defined options and passes arguments on to the individual method scripts.

# TODO: add validity check for DATA_PATH

# Load defaults
source ASsembly_defaults.conf


# Override defaults with given options
while getopts "c:g:hm:" opt; do
case $opt in
	c)
	NCORES=$OPTARG
	;;
	g)
	GTF=$OPTARG
	;;
	m)
	METHODS=$OPTARG
	;;
	o)
	OUTPUT_PATH=$OPTARG
	;;
	h)
	echo "Sorry, help text not available yet!"
	exit 0
	;;
esac
done



DATA_PATH=${@:$OPTIND:1} # First argument should be path to folder storing bam-files folder (and optionally other files)

# If no output path is provided, set data path as output path
if [ -z ${OUTPUT_PATH+x} ]; then OUTPUT_PATH=$DATA_PATH ;  echo "Using default output folder: $OUTPUT_PATH" ; fi

# Add folders to output path structure: int, toGSEA
OUTPUT_PATH=${OUTPUT_PATH}/int
mkdir -p $OUTPUT_PATH
mkdir -p $OUTPUT_PATH/toGSEA

# Create log-file in output path
LOG_FILE="${OUTPUT_PATH}/ASsembly_run_$(date +%Y-%m-%d_%H:%M).log"
cat > $log
exec 3>&1 1>>${LOG_FILE} 2>&1
#cat > $log
#exec 3>&1 1>>${log} 2>&1
#exec > >(tee -i ${log})
#exec 2>&1

# If no gtf-file was supplied, attempt to find gtf-file on data path
if [ -z ${GTF+x} ]; then echo "Retrieving GTF-file from defined data path" | tee /dev/fd/3
nGTF=$(find ${DATA_PATH}/feature_files -name '*.gtf' | wc -l)
	if [ $nGTF != "1" ]; then
		echo " Error: $nGTF gtf-files in feature_files directory. Make sure only 1 is in the dedicated folder or provide a gtf-file manually, using option '-g'." | tee /dev/fd/3
		exit 1
	fi
GTF=$(find ${DATA_PATH}/feature_files -name '*.gtf')
fi
#if [ -z ${+x} ]; then echo "" ; ; fi  # Line for adding more catches of unsupplied parameters


# make comma-separated lists of control and treatment bam-files
cBams=$(find $DATA_PATH/bam_files *.bam | grep "$CONTROL"| tr '\n' ',' | sed 's/.$//')
tBams=$(find $DATA_PATH/bam_files *.bam | grep "$TREATMENT"| tr '\n' ',' | sed 's/.$//')
echo $cBams | tee /dev/fd/3
echo $tBams | tee /dev/fd/3


echo "$(date +%Y-%m-%d_%H:%M): Start individual analyses" | tee /dev/fd/3
IFS=',' read -r -a METHODS <<< "$METHODS"
for methodName in "${METHODS[@]}";do
case $methodName in
	dexseq)
	echo "$(date +%Y-%m-%d_%H:%M): Start DEXSeq_analysis" | tee /dev/fd/3
	bash DEXSeq_master.sh $DATA_PATH $NCORES $GTF $OUTPUT_PATH | tee /dev/fd/3
	echo "$(date +%Y-%m-%d_%H:%M): DEXseq done" | tee /dev/fd/3
	;;
	cufflinks)
	echo "$(date +%Y-%m-%d_%H:%M): Start cufflinks analysis" | tee /dev/fd/3
	bash cufflinks_master.sh $DATA_PATH $NCORES $GTF $OUTPUT_PATH $tBams $cBams | tee /dev/fd/3
	echo "$(date +%Y-%m-%d_%H:%M): Cufflinks done" | tee /dev/fd/3
	;;
	rmats)
	echo "$(date +%Y-%m-%d_%H:%M): Start rMATS analysis" | tee /dev/fd/3
	bash rMATS_master.sh $DATA_PATH $NCORES $GTF $OUTPUT_PATH $tBams $cBams $READLENGTH $PAIRED | tee /dev/fd/3
	echo "$(date +%Y-%m-%d_%H:%M): rMATS done" | tee /dev/fd/3
	;;
	seqgsea) echo "$(date +%Y-%m-%d_%H:%M): Start SeqGSEA analysis" | tee /dev/fd/3
	echo "$(date +%Y-%m-%d_%H:%M): SeqGSEA done" | tee /dev/fd/3
	;;
	aspli) echo "$(date +%Y-%m-%d_%H:%M): Start ASpli analysis" | tee /dev/fd/3
	echo "$(date +%Y-%m-%d_%H:%M): ASpli done" | tee /dev/fd/3
	;;
	*) echo "Error: method $methodName unknown." | tee /dev/fd/3
	;;
esac
done
