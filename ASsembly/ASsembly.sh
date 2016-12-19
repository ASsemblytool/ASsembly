#!/bin/bash

#for debugging purposes
#echo "Press CTRL+C to proceed."
#trap "pkill -f 'sleep 1h'" INT
#trap "set +x ; sleep 1h ; set -x" DEBUG

# ASsembly script
#
# Assures input is valid, applies any defined options and passes arguments on to the individual method scripts.

# Load defaults
source ASsembly_defaults.conf
FIND_THRESHOLDS=false

# Override defaults with given options
while getopts "c:t:g:o:p:Th" opt; do
case $opt in
	c) CONTROL=$OPTARG;;
	t) TREATMENT=$OPTARG;;
	g) GTF=$OPTARG;;
	o) OUTPUT_PATH=$OPTARG;;
	p)
	echo "Loading supplied config file"
	source $OPTARG
	;;
	T)
	echo "Starting threshold determination routine"
	FIND_THRESHOLDS=true
	;;
	h)
	echo "Sorry, help text not available yet!"
	exit 0
	;;
esac
done


######PARSE ARGUMENTS######
# CONTROL (-c), TREATMENT (-t), DATA_PATH (positional arg 1), GTF file (positional arg 2)
DATA_PATH=${@:$OPTIND:1} # Folder storing bam-files folder (and optionally other files)
GTF=${@:$OPTIND+1:1} # path and file name of GTF-file
if [ -z "$DATA_PATH" ]; then echo "Error: supply a data path"; exit 1; fi
if [ -z "$GTF" ]; then echo "Error: supply a gtf-file"; exit 1; fi
#TODO: add validity checks for all arguments



# Folder management
INT=$OUTPUT_PATH/int
GSEA=$INT/GSEA
mkdir -p $INT
if $FIND_THRESHOLDS; then		# To Make sure that thresholds are gone when thresholds are searched
	rm -rf $INT/GSEA
	mkdir $GSEA
	rm -rf $INT/LM
fi

# Create log-file in output path
LOG_FILE="$INT/ASsembly_run_$(date +%Y-%m-%d_%H:%M).log"
exec 3>&1 1>>${LOG_FILE} 2>&1

#if [ -z ${+x} ]; then echo "" ; ; fi  # Line for adding catches of unsupplied parameters

# make comma-separated lists of control and treatment bam-files
cBams=$(find $DATA_PATH/bam_files -name *.bam | grep "$CONTROL"| tr '\n' ',' | sed 's/.$//')
tBams=$(find $DATA_PATH/bam_files -name *.bam | grep "$TREATMENT"| tr '\n' ',' | sed 's/.$//')
allBams=$cBams","$tBams

######INDIVIDUAL ANALYSES######
echo "$(date +%Y-%m-%d_%H:%M): Start individual analyses" | tee /dev/fd/3
METHODS_BACKUP=$METHODS
IFS=',' read -r -a METHODS <<< "$METHODS"
IFS=',' read -r -a TYPES <<< "$TYPES"
for methodName in "${METHODS[@]}";do
case $methodName in
	DEXSeq)
	echo "$(date +%Y-%m-%d_%H:%M): Start DEXSeq_analysis" | tee /dev/fd/3
	bash DEXSeq_master.sh $DATA_PATH $NCORES $GTF $TREATMENT $CONTROL $INT | tee /dev/fd/3
	echo "$(date +%Y-%m-%d_%H:%M): DEXseq done" | tee /dev/fd/3
	;;
	Cufflinks)
	echo "$(date +%Y-%m-%d_%H:%M): Start cufflinks analysis" | tee /dev/fd/3
	bash cufflinks_master.sh $DATA_PATH $NCORES $GTF $INT $tBams $cBams | tee /dev/fd/3
	echo "$(date +%Y-%m-%d_%H:%M): Cufflinks done" | tee /dev/fd/3
	;;
	MATS)
	echo "$(date +%Y-%m-%d_%H:%M): Start rMATS analysis" | tee /dev/fd/3
	bash rMATS_master.sh $DATA_PATH $NCORES $GTF $INT $tBams $cBams $READLENGTH $PAIRED | tee /dev/fd/3
	echo "$(date +%Y-%m-%d_%H:%M): rMATS done" | tee /dev/fd/3
	;;
	SeqGSEA) echo "$(date +%Y-%m-%d_%H:%M): Start SeqGSEA analysis" | tee /dev/fd/3
	echo "$(date +%Y-%m-%d_%H:%M): SeqGSEA done" | tee /dev/fd/3
	;;
	ASpli) echo "$(date +%Y-%m-%d_%H:%M): Start ASpli analysis" | tee /dev/fd/3
	bash ASpli_masterscript.sh $DATA_PATH $GTF $INT $CONTROL $TREATMENT
	echo "$(date +%Y-%m-%d_%H:%M): ASpli done" | tee /dev/fd/3
	;;
	*) echo "Error: method $methodName unknown." | tee /dev/fd/3
	;;
esac
done

######RUN SUPPA######
mkdir -p $INT/SuppaEvents
python3.4 $SUPPA_PATH/suppa.py generateEvents -i $GTF -o $INT/SuppaEvents/ -e SE SS RI MX FL | tee /dev/fd/3
for file in $INT/SuppaEvents/*
        do
                mv "$file" "${file/_/}"
        done

echo "Suppa alternative splicing events lists generated." | tee /dev/fd/3

cat $INT/SuppaEvents/*.ioe | cut -f3 | cut -d ":" -f1 |  awk '{gsub(/;/,"\t");print}' | sort | grep -v "event_id" > $INT/allEventsgenome.txt

######CUTOFF/GSEA PREPARATION######
for methodName in "${METHODS[@]}";do
	case $methodName in
		# For methods that do not autodetect AS type
		Cufflinks|DEXSeq|SeqGSEA)
		        python GSEAconverter_genes.py -f $GSEA/$methodName/*.txt -el $INT/allEventsgenome.txt -m $methodName -d $GSEA/$methodName/
	        ;;
		# For methods that do autodetect AS type
		MATS|ASpli)
			python GSEAconverter_events.py $GSEA/$methodName/*.txt
		;;
	esac
done

######START THRESHOLD DETERMINATION & MODEL TRAINING######
if $FIND_THRESHOLDS; then
######CHECK GENES IN BAM FILES######
	rm -rf $INT/SimulatedGenes
	mkdir -p $INT/SimulatedGenes
	
	for file in $(echo $allBams | sed "s/,/ /g"); do
		fileNoPath=$(echo ${file##*/})
		samtools view -h $file | grep -v "@" | cut -d ":" -f3 | cut -d "." -f1 | sort | uniq > "$INT/SimulatedGenes/"$fileNoPath".txt"
	done
	
	cat $INT/SimulatedGenes/* | sort | uniq > $INT/SimulatedGenes/GenesinBams.txt
	ngenes=$(wc -l $INT/SimulatedGenes/GenesinBams.txt | awk '{print $1}')
	echo 'There are '$ngenes' genes expressed in the bam files' | tee /dev/fd/3
	
	#####GENERATE TRUE GMT FILE####
	echo "$(date +%Y-%m-%d_%H:%M): Generating true events lists" | tee /dev/fd/3
	grep -f $INT/SimulatedGenes/GenesinBams.txt $DATA_PATH/DiffASgenes.txt > $INT/TrueEventlist.txt
	
	mkdir -p $INT/gmts
	
	for type in "${TYPES[@]}";do
		echo ${type}$(grep $type $INT/TrueEventlist.txt | cut -f1 | awk '{for (f=1;f<=NF;f++) col[f] = col[f]":"$f} END {for (f=1;f<=NF;f++) print col[f]}') | tr ':' '	' > $INT/gmts/${type}.gmt
	done
	cat $INT/gmts/*.gmt > $INT/trueEvents.gmt
	
	######GSEA THRESHOLD DETERMINATION######
	echo "$(date +%Y-%m-%d_%H:%M): Determining GSEA thresholds" | tee /dev/fd/3
	echo "Thresholds" > $INT/Thresholds.txt
	
	for methodName in "${METHODS[@]}";do
	Rscript GSEA.R -d $GSEA/$methodName --se SE.rnk --ri RI.rnk --a3 A3.rnk --a5 A5.rnk --mx MX.rnk --al AL.rnk --af AF.rnk --rgmt $INT/trueEvents.gmt -m $methodName
	
	   	echo "$methodName" >> $INT/Thresholds.txt
	        cut -f1,5 $GSEA/$methodName/${methodName}Statistics.txt  >> $INT/Thresholds.txt
	done
	######LINEAR MODEL TRAINING######
	echo "$(date +%Y-%m-%d_%H:%M): Training linear model for results ensemble" | tee /dev/fd/3
	Rscript ASsemblyLM.R construct $INT $METHODS_BACKUP $LM_PRIORITY $LM_THRESHOLD
	
	######TREE MODEL TRAINING######
	
	################################
else
	######THRESHOLD & MODEL APPLICATION######
	for methodName in "${METHODS[@]}";do
		for asType in "${TYPES[@]}"; do
			nLines=($(wc -l ${GSEA}/${methodName}/${asType}.rnk))
			curThres=$(grep $asType ${GSEA}/${methodName}/${methodName}Statistics.txt| cut -f5)
			nTrues=$(awk "BEGIN { pc = $nLines*$curThres; i=int(pc); print (pc-i<0.5)?i:i+1 }")
			tail -q $GSEA/$methodName/${asType}.rnk -n $nTrues > $GSEA/$methodName/LeadingEdge_genes/LE${asType}genes.txt
		done
	done
	Rscript ASsemblyLM.R predict $INT $METHODS_BACKUP $LM_PRIORITY $LM_THRESHOLD
fi

######
