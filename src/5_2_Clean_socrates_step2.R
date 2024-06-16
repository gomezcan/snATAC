###################
## Step 2: post-QC
###################
suppressWarnings(library(Socrates))
suppressWarnings(library(tidyverse)) 
suppressWarnings(library(ggplot2))
suppressWarnings(library(data.table))
suppressWarnings(library(patchwork))

args <- commandArgs(T)


################################
####       Functions        ####
################################

chop=function(myStr,mySep,myField){
  
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

Plot_QC_cell <- function(soc.obj, threshold){
  # estimate log10 number of accessible regions per cell
  # soc.obj = obj_full
  cell.counts <- Matrix::colSums(soc.obj$counts) 
  cell.counts <- tibble(bc=names(cell.counts), count=as.vector(cell.counts))
  
  # estimate peak accessibility frequency across cells
  site.freq <- Matrix::rowMeans(soc.obj$counts)
  site.freq <- tibble(peak=names(site.freq), freq=as.vector(site.freq))
  
  
  cell.counts %>%
    ggplot(aes(x=count))+ 
    geom_density() +
    scale_x_log10(label=comma) +
    annotation_logticks(sides = "b", color = 'black') +
    geom_vline(xintercept = threshold) +
    ylab('Density') + 
    xlab(bquote(log[10]~" bins x cell")) +
    theme_bw() -> Plot_cell_counts
  
  site.freq %>%
    ggplot(aes(x=freq))+ 
    geom_density() +
    #scale_x_log10(label=comma) +
    #annotation_logticks(sides = "b", color = 'black') +
    #geom_vline(xintercept = 1000) +
    xlab('Cell fraction x bin') + 
    ylab("Frequence") +
    theme_bw() -> Plot_peak_freq
  
  Plot_1 <- Plot_cell_counts | Plot_peak_freq
  
  return(Plot_1)
}

################################


# Read rds before soc object created
# Sample = "Sample_P4.Genotype_genome"
sample <- as.character(args[1]) # Sample_P4.Genotype_genome

# Read object 
obj_full <-  readRDS(paste0(sample,".raw.before.soc.rds"))

print(paste0(" Reading .rds with ", nrow(obj_full$meta), " cell .."))
#obj$gff <- "Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Mt.Pt.gff3.gz"
# obj_m$bedpath <- "../7_socrates/qc.Sample_S2.tn5.mq10.bed.gz"
# obj_m$annpath <- "../7_socrates/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.Mt.Pt.gff3.gz"
# obj_m$chrpath <- "../7_socrates/ChrSize"

# generate sparse matrix
obj_full <- generateMatrix(obj_full, filtered=F,  windows=500, peaks=F, verbose=T)

# convert to Socrates format for downstream analysis. 
obj_full <- convertSparseData(obj_full, verbose=T)

print(paste0(" Created sparse matrix for ", sample, " .."))

# Add metadata for cells QCq
obj_full$meta[,"Sample"] <- chop(obj_full$meta$cellID, '[_]',2)
obj_full$meta[,"Genotype"] <- chop(obj_full$meta$cellID, '[_]',3)



## Filter matrix 
### Estimate log10 number of accessible regions per cell
#### Plot before
Plot_QC_before <- Plot_QC_cell(obj_full, 100)
#Plot_QC_before <- Plot_QC_before + plot_annotation( title = 'Before Filter')

#### Clean
obj_full_qc <- cleanData(obj_full, min.c=100, verbose=T)
#obj_full_qc2 <- cleanData(obj_full, min.c=1000, min.t=0.001, max.t=0.005, verbose=T)

##### Filter matrix 
# save object for integration with other libraries
saveRDS(obj_full_qc, file=paste0(sample,'.soc.clean.rds'))
#####

#### Plot after
Plot_QC_after  <- Plot_QC_cell(obj_full_qc, 100)
#Plot_QC_after <- Plot_QC_after + plot_annotation( title = 'After filter: min.c=100')


lable1 <- paste0("Before\nFilter\nTotal.bc: ", nrow(obj_full$meta))
lable2 <- paste0("After filter:\nmin.c=100\nTotal.bc: ", nrow(obj_full_qc$meta))

Plot_QC <- (Plot_QC_before/Plot_QC_after) + 
  plot_annotation( tag_levels = list(c(lable1,"", lable2, "")))

pdf(paste0("Plots/Summary_qc_",sample,".pdf"), width = 6, height = 6)
print(Plot_QC)
dev.off()

print(paste0(".. Done QC plots .."))

# Plot total cell per genotype
Cellg <- obj_full_qc$meta

Cellg <- obj_full_qc$meta$Genotype %>% table %>% as_tibble() %>%
  dplyr::rename(Genotype = '.') %>%
  ggplot(aes(x=reorder(Genotype, -n), y=n)) +
  geom_bar(stat="identity", fill='#FFCCFF')+
  xlab('Genotype assigment')+
  ylab('Cells') +
  labs(subtitle=gsub('Sample_', '', sample)) +
  geom_text(aes(x=reorder(Genotype, -n), y=n, label=n), size=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

pdf(paste0("Plots/FinalCell_", sample, ".pdf"), width = 3, height = 2)
print(Cellg)
dev.off()


print(paste0(".. Done final plots .."))
