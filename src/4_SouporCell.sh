#!/bin/bash

ml Bioinformatics samtools/1.13-fwwss5n vcftools/0.1.15  bcftools/1.12-g4b275e singularity

runDemux(){
	# inputs
	i=$1
	vcf=$2

	# New inputs
	bam=qc.$i.BC.mq10.rmdup.bam
	id=Pool_$i
	out=$id

	# Sample info
	samples=$Pools_DB_list
	fasta=Zm-B73-REFERENCE-NAM-5.0.chrs.mt.pt.fa
		  	
	# Make output dir
	if [ ! -d $out ]; then
		mkdir $out
	fi
	
	# Subset genotypes
	if [ ! -f $out/$id.INPUT.genotypes.txt ]; then
		grep $id $samples | cut -f1 - | uniq - > $out/$id.INPUT.genotypes.txt
	fi
	
	# filter VCF
	if [ ! -f $out/$id.filtered.vcf ]; then
		bcftools view --thread 20 --genotype het -S $out/$id.INPUT.genotypes.txt $vcf -o $out/$id.filtered.vcf;	
	fi

	# subset barcodes
	if [ ! -f $out/$id.BARCODES.txt ]; then
		zcat qc.${i}.tn5.mq10.bed.gz | cut -f4 | sort | uniq -c | awk '{if($1>=100) print $2}' | cut -d':' -f3 > $out/$id.BARCODES.txt
		#cat $i/outs/filtered_peak_bc_matrix/barcodes.tsv | cut -d'-' -f1 >  $out/$id.BARCODES.txt
	fi

	# save genotypes to var
	genos=$(while read i; do
		echo -n "$i "
		done < $out/$id.INPUT.genotypes.txt)
	
	num_k=$(wc -l < $out/$id.INPUT.genotypes.txt)
	echo " - k = $num_k | $genos"

	# 			--skip_remap True \
	# --known_genotypes $out/$id.filtered.vcf #--known_genotypes_sample_names $genos;
	singularity exec ~/souporcell_latest.sif souporcell_pipeline.py --bam $bam \
			--no_umi True \
			--barcodes $out/$id.BARCODES.txt \
			--fasta $fasta --threads 20 \
			--out_dir $out -k $num_k \
			--min_alt 5 --min_ref 5 
	
	# Demuxafy correlation vcf used and and reference snps
	singularity exec ~/Demuxafy.sif Assign_Indiv_by_Geno.R -r $vcf \
		-c $out/cluster_genotypes.vcf -o $out

}

export -f runDemux

# runDemux Pool Name vcf_file
# runDemux $1 $2
# runDemux $1 ../hmp321_agpv5.chr.sorted.reduced_MAF_hets.NAM.vcf.gz

Pools_DB_list=$1
# Pools_DB_list="Pools_DB.txt"

parallel -j 6 runDemux ::: $(cut -f2 $Pools_DB_list | sort -u | sed 's/Pool_//g') ::: hmp321_agpv5.chr.sorted.reduced_MAF_hets.NAM.vcf.gz
