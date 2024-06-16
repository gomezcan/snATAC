## process 282 maize leaf data ##

# load libraries
suppressMessages(library(Socrates))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

# arguments
args <- commandArgs(T)
if(length(args) != 5){stop("Rscript QC_scATAC_data.R <bed> <prefix> <gff/gtf> <chr.fai>")}

# load data
bed <- as.character(args[1])
out <- as.character(args[2])
ann <- as.character(args[3])
chr <- as.character(args[4])

Genotyping <- as.character(args[5]) # type: "chr1" (defaul) or "genome"

# bed <- "qc.Sample_P4.tn5.mq10.bed.gz"
# out <- "Sample_P4"
# ann <- "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Mt.Pt.gff3.gz"
# chr <- "ChrSize"
# Genotyping = 'genome'
#rm(list=ls())
############################
#### Summary Genotyping ####
############################
# load cluster calls from genotyping based on chr1

if (Genotyping=="genome") {
  clusters <- as.data.table(read.table(paste0('Genotyping_Summary_', out,".txt"), sep = '\t', header = T))
  print(' .. using genotyping from whole genome')
} else {
  chr_out <- paste0(out,'_chr1')
  clusters <- as.data.table(read.table(paste0('Genotyping_Summary_', chr_out,".txt"), sep = '\t', header = T))
  print(' .. using genotyping from chr1')
}


############################

############################
#### Filter out doublet ####
singlets <- subset(clusters, status=='singlet')

# rename
singlets$barcode <- paste0('CB:Z:',singlets$barcode)
singlets[,'barcode_id'] <- paste0(singlets$barcode,"_", chop(out, '[_]', 2), '_', singlets$assignment)


# read bed files
bed_file <- read.table(bed, sep = '\t', header = F)
nrow(bed_file)

# remove doublets
#nrow(bed_file)
bed_file <- subset(bed_file, V4 %in% singlets$barcode)
#nrow(bed_file)
bed_file <- left_join(bed_file, singlets[,c("barcode","barcode_id")], by=c('V4'='barcode'))
bed_file$V4  <- bed_file$barcode_id 
bed_file <- bed_file[,1:5]


if (Genotyping=="genome") {
 	write.table(bed_file, paste0("singlets.Genotype_genome.", gsub('.gz', '', bed)), sep = '\t', quote = F, row.names = F, col.names = F)
	system(paste("pigz", paste0("singlets.Genotype_genome.", gsub('.gz', '', bed))))
	
	# redefine bed path
	bed <- paste0("singlets.Genotype_genome.", bed)

} else {
  	write.table(bed_file, paste0("singlets.Genotype_chr.", gsub('.gz', '', bed)), sep = '\t', quote = F, row.names = F, col.names = F)
	system(paste("pigz", paste0("singlets.Genotype_chr.", gsub('.gz', '', bed))))
	
	# redefine bed path
	bed <- paste0("singlets.Genotype_chr.", bed)
}

rm(bed_file)

# load objects
obj <- loadBEDandGenomeData(bed, ann, chr, attribute="transcript_id")

# count organellar reads
obj <- countRemoveOrganelle(obj, org_scaffolds=c("Mt", "Pt"), remove_reads=T)


# call ACRs
if (Genotyping=="genome") {
 	out_dir 	<-  paste0(out,"_genome_peaks")
 	out_dir_tem <- paste0(out, '_genome_macs2_temp')

} else {
 	out_dir 	<-  paste0(out,"_chr_peaks")
 	out_dir_tem <- paste0(out, '_chr_macs2_temp')
}

obj <- callACRs(obj, genomesize=1.6e9, 
                shift= -75, 
                extsize=150,
                fdr=0.1,
                output=out_dir, 
                tempdir=out_dir_tem, 
                verbose=T)

# build metadata
obj <- buildMetaData(obj, tss.window=2000, verbose=TRUE)

# save QC object
if (Genotyping=="genome") { 	
 	saveRDS(obj, file=paste0(out,".Genotype_genome.raw.before.soc.rds"))
 	
} else {
	saveRDS(obj, file=paste0(out,".Genotype_chr.raw.before.soc.rds"))

}



# Generate sparse matrix by bins
obj <- generateMatrix(obj, filtered=F,  windows=500, peaks=F, verbose=T)

# add metadata for cells QC
#obj <- isCell(obj)
