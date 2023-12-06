#!/bin/bash

###################################
#######       Modules     #########
###################################

ml Bioinformatics cellrangerATAC/2.1.0

###################################
###    Command lines to Run     ###
###################################

Samples=$1 # Sample id file
Index=$(Index_genome_CellrangerATAC)

while read line; do 
	#
	cellranger-atac count --id=Sample_${line} --reference=$Index --fastqs ${line} --sample $line;
done < $Samples
