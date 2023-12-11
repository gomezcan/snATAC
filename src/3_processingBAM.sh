#!/bin/bash

###################################
#######       Modules     #########
###################################
ml Bioinformatics samtools/1.13-fwwss5n picard-tools/3.0.0

###################################
###    Command lines to Run     ###
###################################

# functions
doCall(){

	# input
	base=$1
	qual=$2
	threads=$3
	
	# filter for mapped reads only
	echo "retaining only mapped reads ..."
	samtools view -@ ${threads} -bhq ${qual} -f 3 ${base}.bam > qc.${base}.mq${qual}.bam

 	# clean BC names
	echo "cleaning BC names - $base ..."
	samtools view -@ 5 -h qc.${base}.mq${qual}.bam | \
		sed 's@\(CB:Z:[ACGT]\{16\}\)-1@\1@' | \
		samtools view -@ 5 -b -o qc.${base}.mq${qual}.bcClean.bam
  
	# run picard
	echo "removing dups - $base ..."
	java -Xmx100g -jar $PICARDLIB/picard.jar MarkDuplicates \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		REMOVE_DUPLICATES=true \
		METRICS_FILE=${base}.metrics \
		I=qc.${base}.mq${qual}.bcClean.bam \
		O=qc.${base}.mq${qual}.rmdup.bam \
		BARCODE_TAG=CB \
		ASSUME_SORT_ORDER=coordinate \
		USE_JDK_DEFLATER=true \
		USE_JDK_INFLATER=true
	
	# fix barcodes
	echo "fixing barcodes and removing multi-mapped reads ..."
	perl fixBC.pl qc.${base}.mq${qual}.rmdup.bam ${base} | samtools view -bhS - > qc.${base}.BC.mq${qual}.rmdup.bam

	# make Tn5 bed files
	echo "making Tn5 bed files ..."
	perl makeTn5bed.pl wc.${base}.BC.mq${qual}.rmdup.bam | sort -k1,1 -k2,2n - | uniq - > qc.${base}.tn5.mq${qual}.bed
	pigz ${base}.tn5.mq${qual}.bed
	

}

export -f doCall

# run processing (change sample1 to the file name without the .bam suffix)
#doCall sample $qual $threads

# Call downloading function
parallel -j 8 doCall :::: SampleList ::: 10 ::: 60
