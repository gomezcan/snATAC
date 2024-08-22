# load libraries
suppressMessages(library(Socrates))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

library(GenomicFeatures)
library(GenomicRanges)
library(Matrix)

rm(list=ls())

###################################
# load arguments
args <- commandArgs(trailingOnly=T)

estGeneActivity <- function(obj, FeatureName="gene", pRange=50000, dRange=100, con=10000){
  
  ## Load Tn5 file with GRanges
  Tn5_Grange <-  GRanges(seqnames = obj$bed$V1,
                         ranges = IRanges(start = obj$bed$V2,
                                          end = obj$bed$V3,
                                          names = obj$bed$V4))
  
  ## Select Features
  if(FeatureName == "gene"){
    Ann <- genes(obj$gff)
  }else if(FeatureName == "transcript"){
    Ann <- transcripts(obj$gff)
  }else{
    message("ERROR: Feature name should be 'gene' or 'transcript")
  }
  
  ## Find overlap in genes
  message(" - building sparse gene body matrix ...")
  hits_Within <- suppressWarnings(findOverlaps(Tn5_Grange,  Ann, minoverlap=1,
                                               type=c("within"), select="all", 
                                               ignore.strand = TRUE))
  
  Intersect <- paste(names(Ann)[hits_Within@to],
                     names(Tn5_Grange)[hits_Within@from],sep="/_Com_/")
  Intersect <- table(Intersect)
  Intersect <- data.frame(geneID = as.factor(as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                                          split="/_Com_/"), "[", 1))),
                          cellID = as.factor(as.character(lapply(strsplit(as.character(rownames(Intersect)),
                                                                          split="/_Com_/"), "[", 2))),
                          activity = as.numeric(as.character(Intersect)))
  obj$gba <- sparseMatrix(i=as.numeric(Intersect$geneID),
                          j=as.numeric(Intersect$cellID),
                          x=Intersect$activity,
                          dimnames=list(levels(Intersect$geneID),
                                        levels(Intersect$cellID)))
  rm(Intersect)
  rm(hits_Within)
  all.genes <- rownames(obj$gba)
  
  # find overlaps in peaks
  message(" - building sparse ACR matrix ...")
  acrs <- GRanges(seqnames=obj$acr$V1,
                  ranges=IRanges(start=obj$acr$V2,
                                 end=obj$acr$V3,
                                 names=paste(obj$acr$V1,obj$acr$V2,obj$acr$V3, sep="_")))
  peak_hits <- suppressWarnings(findOverlaps(Tn5_Grange, acrs, 
                                             minoverlap=1, 
                                             type=c("within"), 
                                             select="all", 
                                             ignore.strand=T))
  
  peak_int <- paste(names(acrs)[peak_hits@to],
                    names(Tn5_Grange)[peak_hits@from],sep="/_Com_/")
  peak_int <- table(peak_int)
  peak_int <- data.frame(acrID = as.factor(as.character(lapply(strsplit(as.character(rownames(peak_int)),
                                                                        split="/_Com_/"), "[", 1))),
                         cellID = as.factor(as.character(lapply(strsplit(as.character(rownames(peak_int)),
                                                                         split="/_Com_/"), "[", 2))),
                         activity = as.numeric(as.character(peak_int)))
  obj$pa <- sparseMatrix(i=as.numeric(peak_int$acrID),
                         j=as.numeric(peak_int$cellID),
                         x=peak_int$activity,
                         dimnames=list(levels(peak_int$acrID),
                                       levels(peak_int$cellID)))
  #message("   peak matrix = ", nrow(obj$pa), " | ", ncol(obj$pa))
  all.peaks <- rownames(obj$pa)
  
  # find peaks within the gene range
  message(" - building ACR x gene body weight matrix ...")
  upstream <- promoters(Ann, upstream=pRange, downstream=0)
  upstream_peaks <- suppressWarnings(findOverlaps(acrs, upstream, minoverlap=1, 
                                                  type=c("within"),
                                                  select="all",
                                                  ignore.strand=T))
  peak_genes <- data.frame(geneID=names(Ann)[upstream_peaks@to],
                           acrID=names(acrs)[upstream_peaks@from])
  all.genes <- unique(c(all.genes), as.character(peak_genes$geneID))
  all.peaks <- unique(c(all.peaks), as.character(peak_genes$acrID))
  peak_genes$distance <- distance(Ann[peak_genes$geneID,],
                                  acrs[peak_genes$acrID,])
  peak_genes$geneID <- factor(as.character(peak_genes$geneID), levels=all.genes)
  peak_genes$acrID <- factor(as.character(peak_genes$acrID), levels=all.peaks)
  peak_genes$weight <- 2^-(peak_genes$distance/con)
  peak_genes <- peak_genes[complete.cases(peak_genes),]
  obj$pgw <- sparseMatrix(i=as.numeric(peak_genes$geneID),
                          j=as.numeric(peak_genes$acrID),
                          x=peak_genes$weight,
                          dims=c(length(levels(peak_genes$geneID)),
                                 length(levels(peak_genes$acrID))),
                          dimnames=list(levels(peak_genes$geneID),
                                        levels(peak_genes$acrID)))
  
  shared <- intersect(rownames(obj$pa), colnames(obj$pgw))
  #message("   peak gene weighted matrix = ", nrow(obj$pgw), " | ", ncol(obj$pgw))
  #message("   number of shared peaks = ", length(shared))
  obj$peak_gene_score <- obj$pgw[,shared] %*% obj$pa[shared,]
  #message("   regulatory activity = ", nrow(obj$peak_gene_score), " | ", ncol(obj$peak_gene_score))
  #message(" - gene body accessibility")
  #print(head(obj$gba[,1:5]))
  #message(" - peak accessibility")
  #print(head(obj$pa[,1:5]))
  #message(" - regulatory activity")
  #print(head(obj$peak_gene_score[,1:5]))
  
  # aggregate regulatory score with gene body score
  shared.cells <- intersect(colnames(obj$peak_gene_score),
                            colnames(obj$gba))
  shared.genes <- intersect(rownames(obj$peak_gene_score),
                            rownames(obj$gba))
  obj$gene_activity <- obj$peak_gene_score[shared.genes,shared.cells] + obj$gba[shared.genes,shared.cells]
  
  return(obj)
}

###################################

#args
obj_name    <- as.character(args[1]) # rds object
bed_file    <- as.character(args[2]) # bed with taq
SampleName  <- as.character(args[3]) # sample name

###################################
#####     load soc obj   ##########
###################################

obj <- readRDS(obj_name)

# add acrs to obj
obj$acr <- read.table("All_Good_BCbyClusters_peaks.08_21_24.narrowPeak")
 
# add bed to obj
#obj$bed <- read.table(gzfile("All.tn5.good.bc.bed.gz"))
obj$bed <- read.table(bed_file)

# add annotation
ann <- "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Mt.Pt.gff3.gz"
gff <- makeTxDbFromGFF(ann, format="gff3", dbxrefTag="Parent")
obj$gff <- gff

# selected top 1000 cell by cluster
head(obj$h.Clusters2)
meta_top <- obj$h.Clusters2 %>%
  dplyr::group_by(LouvainClusters) %>%
  dplyr::top_n(500, log10nSites)

# top cells id
bc_to_select <- meta_top$cellID

#
obj$h.Clusters2 <- obj$h.Clusters2[obj$h.Clusters2$cellID %in% bc_to_select,]

# filter bed
obj$bed <- obj$bed[as.character(obj$bed$V4) %in% bc_to_select,]


## Summary plots 
obj <- estGeneActivity(obj, FeatureName="gene", pRange=500, con=5000)


# save data --------------------------------------------------------------
saveRDS(obj, file=paste0("GeneActivity/Full_obj.", SampleName, '.rds'))
saveRDS(obj$gene_activity, file=paste0("GeneActivity/Gene_activity.", SampleName, '.rds'))
saveRDS(obj$peak_gene_score, file=paste0("GeneActivity/Peak_gene_score.", SampleName, '.rds'))

saveRDS(obj$gba, file=paste0("GeneActivity/bga_/gba.", SampleName, '.rds'))
saveRDS(obj$pa, file=paste0("GeneActivity/pa_/pa.", SampleName, '.rds'))
saveRDS(obj$pgw, file=paste0("GeneActivity/pgw_/pgw.", SampleName, '.rds'))
