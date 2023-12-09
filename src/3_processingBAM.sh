#!/bin/bash

###################################
#######       Modules     #########
###################################
ml Bioinformatics samtools/1.13-fwwss5n picard-tools/3.0.0

###################################
###    Command lines to Run     ###
###################################

# threads
threads=5
qual=10

# functions
doCall(){
	# input paramters 
	base=$1
	qual=$2
	threads=$3

	# filter for mapped reads only
	echo "retaining only mapped reads ..."
	samtools view -@ $threads -bhq $qual -f 3 $base.bam > $base.mq$qual.bam

	# run picard
	echo "removing dups - $base ..."
	java -Xmx100g -jar $PICARDLIB/picard.jar MarkDuplicates \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		REMOVE_DUPLICATES=true \
		METRICS_FILE=$base.metrics \
		I=$base.mq$qual.bam \
		O=$base.mq$qual.rmdup.bam \
		BARCODE_TAG=CB \
		ASSUME_SORT_ORDER=coordinate \
		USE_JDK_DEFLATER=true \
		USE_JDK_INFLATER=true

	# fix barcodes
	echo "fixing barcodes and removing multi-mapped reads ..."
	perl fixBC.pl $base.mq$qual.rmdup.bam $base | samtools view -bhS - > $base.BC.mq$qual.rmdup.bam

	# make Tn5 bed files
	echo "making Tn5 bed files ..."
	perl makeTn5bed.pl $base.BC.mq$qual.rmdup.bam | sort -k1,1 -k2,2n - | uniq - > $base.tn5.mq$qual.bed
	gzip $base.tn5.mq$qual.bed

}
export -f doCall

# run processing (change sample1 to the file name without the .bam suffix)
doCall sample1 $qual $threads
