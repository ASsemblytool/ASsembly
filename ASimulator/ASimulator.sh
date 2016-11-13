#!/bin/bash

# ASimulator (v1.0): Simulation of Differential Alternative Splicing 
# Authors: Carmen Bravo, Carlos de Lannoy, Emma Houben
# Email: assemblytoolkuleuven@gmail.com
# Github: https://github.com/ASsemblytool/ASsembly

# 0. Help function (usage)
programname=$0

function usage {
    echo "usage: $programname [-hn] [-e simpleEvents | complexEvents | allEvents] [--gtf] [-t numberOfASgenes] [-r numberOfReplicates] [--palt1 paltRatio] [--palt2 paltRatio] [-m meanbasecoverage] [--gsnap fullcDNAlibrary] --group1 identifierGroup1 --group2 identifierGroup2"
    echo "Authors: Carmen Bravo, Carlos de Lannoy, Emma Houben."
    echo "Email: assemblytoolkuleuven@gmail.com"
    echo "Github: https://github.com/ASsemblytool/ASsembly"
    echo "Version: 1.0"
    echo "  -e simpleEvents | complexEvents | allEvents			specify which genes will be differentially spliced according to event type. Default is random selection."
    echo "  --group1 identifierGroup1					(mandatory) specify identifier of group 1."
    echo "  --group2 identifierGroup2					(mandatory) specify identifier of group 2."
    echo "  --gsnap fullcDNAlibrary 					turn on gsnap alignment (cDNA library file must be provided as input). Default is false."
    echo "  --gtf								create gtf annotation. Default is false."
    echo "  -h								display help."
    echo "  -m meanbasecoverage						specify mean base coverage. Default is 25."
    echo "  -n								produce an incomplete annotation file for novel discovery testing. Default is false."
    echo "  --palt1 paltRatio						specify PALT ratio of group 1. Default is 0.2."
    echo "  --palt2 paltRatio						specify PALT ratio of group 2. Default is 0.6."
    echo "  -r numberOfReplicates						specify number of biological replicates. Default is 3."
    echo "  -t numberOfASgenes						specify number of genes to be differentially spliced between groups if random selection is chosen. Default is 2000."
    exit 1
}

# 1. Define default arguments
gtfflag='false'
eventTypeflag=''
novelflag='false'
fgpalt=0.2
sgpalt=0.6
ntargflag=2000
replflag=3
gsnapflag=''
meanbasecoverage=25
fgflag=''
fgflag=''

# 2. Parse arguments
# NOTE: MAC getopts does not accept long options. The solution is to convert them to 
# single character options.
for arg in "$@"; do
  shift
  case "$arg" in
    "--gtf") set -- "$@" "-g" ;;
    "--gsnap") set -- "$@" "-a" ;;
    "--group1") set -- "$@" "-1" ;;
    "--group2") set -- "$@" "-2" ;;
    "--palt1") set -- "$@" "-f" ;;
    "--palt2") set -- "$@" "-s" ;;
    *)        set -- "$@" "$arg"
  esac
done

OPTIND=1 
while getopts 'e:t:1::2::r:f:s:m:a:hgn' flag; do
  case "${flag}" in
    e) eventTypeflag="${OPTARG}" ;;
    t) ntargflag="${OPTARG}" ;;
    r) replflag="${OPTARG}" ;;
    1) fgflag="${OPTARG}" ;;
    2) sgflag="${OPTARG}" ;;
    f) fgpalt="${OPTARG}" ;;
    s) sgpalt="${OPTARG}" ;;
    m) meanbasecoverage="${OPTARG}" ;;
    a) gsnapflag="${OPTARG}" ;;
    n) novelflag="true" ;;
    g) gtfflag="true" ;;
    h) usage
  esac
done
shift $(expr $OPTIND - 1) # Not really necessary since no argument is taken by position.

# 3. Check that group identifiers have been given.
if [ "$fgflag" == '' ] || [ "$sgflag" == '' ]
	then
		echo 'Error: Provide group identifiers.'
		exit
fi

# 4. Define files paths (retrieved from myPaths.path)
MASTERDIR=$(grep 'MASTERDIR' myPaths.path | cut -d "=" -f2)
FUNCTIONS=$(grep 'FUNCTIONS' myPaths.path | cut -d "=" -f2)
INT=$(grep 'INT' myPaths.path | cut -d "=" -f2)
ANNOT=$(grep 'ANNOT' myPaths.path | cut -d "=" -f2)
GFF_ANN=$(grep 'GFF_ANN' myPaths.path | cut -d "=" -f2)
GTF_ANN=$(grep 'GTF_ANN' myPaths.path | cut -d "=" -f2)
BAM_FILES=$(grep 'BAM_FILES' myPaths.path | cut -d "=" -f2)
FILT_BAMFILES=$(grep 'FILT_BAMFILES' myPaths.path | cut -d "=" -f2)
REF_FASTA=$(grep 'REF_FASTA' myPaths.path | cut -d "=" -f2)
FSPARA=$(grep 'FSPARA' myPaths.path | cut -d "=" -f2)
SIMBAM=$(grep 'SIMBAM' myPaths.path | cut -d "=" -f2)

# 5. Run AS finder
python $FUNCTIONS/ASfinder.py $GFF_ANN -o $INT/1-ASgenomeReg.bed
echo "Regions corresponding to alternatively spliced genes found."

# 6. Filter bam files and define control/treatment group
if [ -d "$FILT_BAMFILES" ]; 
	then 
  		rm -r "$FILT_BAMFILES" 
fi
mkdir "$FILT_BAMFILES"
for file in $BAM_FILES/*.bam
 	do
		echo $file
		samtools view -b -f 0x2 -L $INT/1-ASgenomeReg.bed $file | samtools view -h | grep -v '\t\*\t' | samtools view -bS - > ${file/$BAM_FILES/$FILT_BAMFILES}
	done  
echo "Bam files have been filtered according to alternatively spliced genes regions."
# 
# 6.a First group files
fg=$(ls $FILT_BAMFILES | grep "$fgflag" | awk '$0="'$FILT_BAMFILES'/"$0' | tr "\n" " ") 
# 6.b Second group files
sg=$(ls $FILT_BAMFILES | grep "$sgflag" | awk '$0="'$FILT_BAMFILES'/"$0' | tr "\n" " ")

echo "The two groups have been defined (group1: "$fgflag" // group2: "$sgflag")." 

# 7. Convert gff to gtf
if [ "$gtfflag" == "true" ]
	then
		gffread $GFF_ANN -T -o $GTF_ANN
		echo "GTF file from GFF annotation has been created."
else
	echo "GTF annotation is provided by user."
fi


8. Generate event list if specified and simulate biological replicates
if [ "$eventTypeflag" == "" ]
	then
		cat $INT/1-ASgenomeReg.bed | cut -f4 | sort | uniq > $INT/3-allASgenesGenome.txt
		echo "No event type simulation specified, proceed with random selection"
		if [ -d "$INT/4-CalNBCounts" ]; 
			then
				rm -r "$INT/4-CalNBCounts"
		fi
		mkdir "$INT/4-CalNBCounts" 
		python $FUNCTIONS/cal_NB_counts.py $GFF_ANN -g1 $fg -g2 $sg -l $ntargflag -n $replflag -o $INT/4-CalNBCounts -of $MASTERDIR/DiffASgenes.txt
		
elif [ "$eventTypeflag" == "simpleEvents" ]
	then
		if [ -d "$INT/2-SuppaEvents" ]; 
			then 
  				rm -r "$INT/2-SuppaEvents" 
		fi
		mkdir $INT/2-SuppaEvents
		python3.4 $FUNCTIONS/suppa/suppa.py generateEvents -i $GTF_ANN -o $INT/2-SuppaEvents/ -e SE SS RI MX FL
		for file in $INT/2-Suppaevents/*
			do 
				mv "$file" "${file/_/}"
			done
		echo "Suppa alternative splicing events lists generated."
		python $FUNCTIONS/simpleEvents.py $INT/2-SuppaEvents/*.ioe -gff $GFF_ANN -o $INT/3-simpleEventsGenome.txt
		sort $INT/3-simpleEventsGenome.txt -o $INT/3-simpleEventsGenome.txt
		echo "Simple events list generated."
		if [ -d "$INT/4-CalNBCounts" ]; 
			then
				rm -r "$INT/4-CalNBCounts"
		fi
		mkdir "$INT/4-CalNBCounts" 
		python $FUNCTIONS/cal_NB_counts.py $GFF_ANN -g1 $fg -g2 $sg -as $INT/3-simpleEventsGenome.txt -n $replflag -o $INT/4-CalNBCounts -of $MASTERDIR/DiffASgenes.txt
		
elif [ "$eventTypeflag" == "complexEvents" ]
	then
		if [ -d "$INT/2-SuppaEvents" ]; 
			then 
  				rm -r "$INT/2-SuppaEvents" 
		fi
		mkdir $INT/2-SuppaEvents
		python3.4 $FUNCTIONS/suppa/suppa.py generateEvents -i $GTF_ANN -o $INT/2-SuppaEvents/ -e SE SS RI MX FL
		for file in $INT/2-Suppaevents/*
			do 
				mv "$file" "${file/_/}"
			done
		echo "Suppa alternative splicing events lists generated."
		python $FUNCTIONS/complexEvents.py $INT/2-SuppaEvents/*.ioe -gff $GFF_ANN -o $INT/3-complexEventsGenome.txt
		sort $INT/3-complexEventsGenome.txt -o $INT/3-complexEventsGenome.txt
		echo "Complex events list generated."
				if [ -d "$INT/4-CalNBCounts" ]; 
			then
				rm -r "$INT/4-CalNBCounts"
		fi
		mkdir "$INT/4-CalNBCounts" 
		python $FUNCTIONS/cal_NB_counts.py $GFF_ANN -g1 $fg -g2 $sg -as $INT/3-complexEventsGenome.txt -n $replflag -o $INT/4-CalNBCounts -of $MASTERDIR/DiffASgenes.txt
		
elif [ "$eventTypeflag" == "allEvents" ]  
	then
		if [ -d "$INT/2-SuppaEvents" ]; 
			then 
  				rm -r "$INT/2-SuppaEvents" 
		fi
		mkdir $INT/2-SuppaEvents
		python3.4 $FUNCTIONS/suppa/suppa.py generateEvents -i $GTF_ANN -o $INT/2-SuppaEvents/ -e SE SS RI MX FL
		for file in $INT/2-Suppaevents/*
			do 
				mv "$file" "${file/_/}"
			done
		echo "Suppa alternative splicing events lists generated."
		python $FUNCTIONS/allEvents.py $INT/2-SuppaEvents/*.ioe -gff $GFF_ANN -o $INT/3-allEventsGenome.txt
		sort $INT/3-allEventsGenome.txt -o $INT/3-allEventsGenome.txt
		echo "Event type list generated."
				if [ -d "$INT/4-CalNBCounts" ]; 
			then
				rm -r "$INT/4-CalNBCounts"
		fi
		mkdir "$INT/4-CalNBCounts" 
		python $FUNCTIONS/cal_NB_counts.py $GFF_ANN -g1 $fg -g2 $sg -as $INT/3-allEventsGenome.txt -n $replflag -o $INT/4-CalNBCounts -of $MASTERDIR/DiffASgenes.txt
else
	echo "Error: Provide a correct event type selection: -e [simpleEvents | complexEvents | allEvents], or default (random selection)."
	exit
fi

echo "Simulation of biological replicates completed."



# 9. Run Flux simulator
export PATH=$PATH:$FUNCTIONS/flux-simulator-1.2.1/bin
cp $FSPARA $FUNCTIONS/flux-simulator-1.2.1/fluxPara.par

if [ "$eventTypeflag" == "" ]
	then
		python $FUNCTIONS/generate_rnaseq.py $INT/4-CalNBCounts/group1.nbcounts $MASTERDIR/DiffASgenes.txt $FUNCTIONS/flux-simulator-1.2.1/fluxPara.par $INT/5-$fgflag/$fgflag -p $fgpalt -c $meanbasecoverage -r $replflag
		python $FUNCTIONS/generate_rnaseq.py $INT/4-CalNBCounts/group2.nbcounts $MASTERDIR/DiffASgenes.txt $FUNCTIONS/flux-simulator-1.2.1/fluxPara.par $INT/6-$sgflag/$sgflag -p $sgpalt -c $meanbasecoverage -r $replflag
		echo "Simulation of differential alternative splicing finished. Check simulated files in $INT/5-$fgflag/$fgflag and $INT/6-$sgflag/$sgflag ."		
else
	cat $MASTERDIR/DiffASgenes.txt | cut -f1 > $INT/4a-DiffASgenesIDs.txt
	if [ -d "$INT/5-$fgflag" ]; 
		then 
   				rm -r "$INT/5-$fgflag" 
 	fi
 	mkdir "$INT/5-$fgflag"
	python $FUNCTIONS/generate_rnaseq.py $INT/4-CalNBCounts/group1.nbcounts $INT/4a-DiffASgenesIDs.txt $FUNCTIONS/flux-simulator-1.2.1/fluxPara.par $INT/5-$fgflag/$fgflag -p $fgpalt -c $meanbasecoverage -r $replflag
	if [ -d "$INT/6-$sgflag" ]; 
		then 
   				rm -r "$INT/6-$sgflag" 
 	fi
 	mkdir "$INT/6-$sgflag"
	python $FUNCTIONS/generate_rnaseq.py $INT/4-CalNBCounts/group2.nbcounts $INT/4a-DiffASgenesIDs.txt $FUNCTIONS/flux-simulator-1.2.1/fluxPara.par $INT/6-$sgflag/$sgflag -p $sgpalt -c $meanbasecoverage -r $replflag
	echo "Simulation of differential alternative splicing finished."
fi

# 10. Generate incomplete annotation if novel isoform discovery is specified

if [ "$novelflag" == 'true' ]
	then 
		python $FUNCTIONS/incompleteAnnotation.py $INT/4a-DiffASgenesIDs.txt -gff $GFF_ANN -o $ANNOT/IncompleteGFF.gff
		gffread $ANNOT/IncompleteGFF.gff -T -o $ANNOT/IncompleteGTF.gtf
		echo "Incomplete annotations generated at $ANNOT"
fi

# 11. GMAP/GSNAP alignment if specified
if [ "$gsnapflag" != '' ]
	then	
		# echo "Starting alignment with gmap/gsnap."
 		gmap_build -d 7-Gmap_genome -D $INT $REF_FASTA/*.fa
  		gmap -d $INT/7-Gmap_genome -f 6 $gsnapflag $INT > $INT/8-splicesites
  		cat $INT/8-splicesites | iit_store -o $INT/9-splicesites_map
		cp $INT/9-splicesites_map.iit $INT/7-Gmap_genome/7-Gmap_genome.maps
 		for i in `seq 1 3`; do
         	file1="$INT/5-$fgflag/"$fgflag"_"$i"_1.fq"
     		file2="${file1%_1.fq}_2.fq"
     		file3="$INT/6-$sgflag/"$sgflag"_"$i"_1.fq"
     		file4="${file3%_1.fq}_2.fq"
     		gsnap -d 7-Gmap_genome -D $INT -s 9-splicesites_map.iit $file1 $file2 -m 10 --format=sam | samtools view -bS > "$SIMBAM/$fgflag"_"$i".bam
     		gsnap -d $INT/7-Gmap_genome -D $INT -s 9-splicesites_map $file3 $file4 -m 10 --format=sam | samtools view -bS > "$SIMBAM/$sgflag"_"$i".bam
 		done
fi
 		
