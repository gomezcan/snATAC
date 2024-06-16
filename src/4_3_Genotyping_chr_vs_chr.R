suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
library(viridis)
library(patchwork)
library(cvms)
library(broom)
library(ggimage)
library("ggnewscale")
library("rsvg")

rm(list=ls())
args <- commandArgs(T)
out <- as.character(args[1])

#############################################
############       Functions     ############
#############################################

chop <- function(myStr,mySep,myField){
  #
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  #
  if(length(myField)>1){
    
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
    
  }
  return(choppedString)
}

Get_GenotypeCor <- function(id_list){
  # made list of cor table by chrs
  # id_list <- CorTables
  corlist <- paste0(id_list, "/ref_clust_pearson_correlations.tsv")
  
  # check if cor matrix available
  corlist <- corlist[unlist(lapply(corlist, file.exists))]
  
  # update list of chrs
  id_list <- chop(corlist, '[/]', 1)
  Sampleid <- unique(chop(id_list, '[_]', 3))
  
  # read cor matrix
  dfcor <- lapply(corlist, fread)
  names(dfcor) <- paste0(chop(id_list, '[_]', 3), ".", chop(id_list, '[_]', 4))
  
  dfcor <- lapply(dfcor, function(x) gather(x, Genotype, Cor, -Cluster))
  dfcor <- rbindlist(dfcor, idcol = T)
  dfcor$Genotype <- gsub('G282set_','', dfcor$Genotype)
  
  return(dfcor)
}

Get_SingletsCounts <- function(id_list){
  # made list of cor table by chrs
  # id_list <- CorTables
  Summarylist <- gsub("Pool", "Genotyping_Summary", id_list)
  Summarylist <- paste0(Summarylist, ".txt")
  # check if cor matrix available
  Summarylist <- Summarylist[unlist(lapply(Summarylist, file.exists))]
  
  # read bc genotyped
  DFSummary <- lapply(Summarylist, function(x) fread(x, select = c("status","assignment")))
  names(DFSummary) <- paste0(chop(Summarylist, '[_]', 4),"_", gsub(".txt", "", chop(Summarylist, '[_]', 5)))
  
  # keep singlets
  DFSummary <- lapply(DFSummary, function(x) subset(x, status =='singlet'))
  
  # Combined DFs
  DFSummary <- rbindlist(DFSummary, idcol = T)
  
  # Counts genotype call call per chr
  DFSummary <- table(DFSummary[,c(".id", "assignment")]) %>%
    as.data.table()
  
  return(DFSummary)
}

Get_Cellgenotypes <- function(id_list){
  # made list of cor table by chrs
  # id_list <- CorTables
  Summarylist <- gsub("Pool", "Genotyping_Summary", id_list)
  Summarylist <- paste0(Summarylist, ".txt")
  #
  # check if cor matrix available
  Summarylist <- Summarylist[unlist(lapply(Summarylist, file.exists))]
  
  # read bc genotyped
  DFSummary <- lapply(Summarylist, function(x) fread(x, select = c("barcode","status","assignment")))
  names(DFSummary) <- paste0(chop(Summarylist, '[_]', 4),"_", gsub(".txt", "", chop(Summarylist, '[_]', 5)))
  names(DFSummary) <- gsub(".txt_NA", "", names(DFSummary))
  
  # Combined DFs
  # keep singlets as only assignment valid
  DFSummary <- rbindlist(DFSummary,idcol = T)
  DFSummary$assignment[DFSummary$status!="singlet"] <- 'unassigned'
  
  return(DFSummary)
  
}


Confusiont_plot <- function(target, prediction){
  
  d_multi <- tibble("target" = CellgenotypesDB[[target]],
                    "prediction" = CellgenotypesDB[[prediction]])
  
  conf_mat <- confusion_matrix(targets = d_multi$target, 
                               predictions = d_multi$prediction)
  
  # plot_confusion_matrix(conf_mat$`Confusion Matrix`[[1]], add_sums = TRUE)
  Plot_ch <- plot_confusion_matrix( 
    conf_mat$`Confusion Matrix`[[1]],
    add_sums = TRUE,
    sums_settings = sum_tile_settings(
      palette = "Oranges",
      label = "Total",
      tc_tile_border_color = "black"),
    font_row_percentages= font(1),
    font_col_percentages= font(1),
    font_normalized=font(2),
    font_counts =font(1.5)
  ) +
  labs(title = prediction) 
  return(Plot_ch)
}
                      
#############################################
                      
#########################################
#### Summary Genotyping:  first part ####
#########################################
# load cluster calls from genotyping
# out <- "Sample_P4"
clusters <- as.data.table(read.table(paste0('Pool_', out, "/clusters.tsv"), sep = '\t', header = T))


# load cluster assignation
clusters_id <-read.table(paste0('Pool_', out, "/ref_clust_pearson_correlations.tsv"), sep = '\t', header = T)
row.names(clusters_id) <- clusters_id$Cluster
clusters_id <- clusters_id[,-c(1)]

# reshape and select highest correlation by genotype
clusters_id <- reshape2::melt(as.matrix(clusters_id)) %>%
  dplyr::group_by(Var1) %>%
  dplyr::top_n(1, value) 


# re-name genotypes
clusters_id$Var2 <- gsub('G282set_', '', clusters_id$Var2)
clusters_id$Var1 <- as.character(clusters_id$Var1)

# replace class
for (j in seq(1:length(clusters$assignment))){
  
  # for is row in genotype, look for each each genotype and replace with name
  list_tem <- data.frame(g=as.character(unlist(strsplit(clusters$assignment[j], '[/]'))))
  
  list_tem <- left_join(list_tem, clusters_id[,1:2], by=c("g"="Var1"))
  
  if (length(list_tem$Var2)>1){
    list_tem <- paste(list_tem$Var2, collapse="/")
  }
  else {
    list_tem <- list_tem$Var2
  }
  clusters$assignment[j] <- list_tem
}

#clusters
print(".. 1. Done genotype name replacement ..")

# Defined high quality sinlges
print(table(clusters$status))

## Reassigment "unassigned" cell based on: 
### 1. log_prob_singleton <= to %75 percentil singlet calls 
singlet_threshold <- quantile(subset(clusters, status=="singlet")$log_prob_singleton, 0.75)

### 1. log_prob_doublet > to %75 percentil singlet calls
doublet_threshold <- quantile(subset(clusters, status=="singlet")$log_prob_doublet, 0.75)

Newsinglet <- subset(clusters, status=="unassigned" & log_prob_singleton <= singlet_threshold & log_prob_singleton > doublet_threshold)

# remove  genotypes unassigned cell
Newsinglet <- Newsinglet[grepl("/", Newsinglet$assignment) == F,]
print(".. 2. Done New singles ..")

# redefine BC in new singlets
clusters$status[clusters$barcode %in%  Newsinglet$barcode] <- 'singlet'
#
write.table(clusters, paste0("Genotyping_Summary_", out, ".txt"), row.names = F, quote = F, sep = '\t')

clusters_Freq <- rbind(
  clusters %>% 
    dplyr::filter(status=='singlet')%>%
    dplyr::select(assignment) %>% table %>% as.data.table() %>% 
    dplyr::mutate(status = 'singlet') %>% dplyr::arrange(-N),
  clusters %>% 
    dplyr::filter(status!='singlet')%>%
    dplyr::select(assignment) %>% table %>% as.data.table() %>% 
    dplyr::mutate(status = 'doublet') %>% dplyr::arrange(-N)
  
)

print(".. 3. Done freq counting ..")
# Plot 1
clusters_Freq %>% 
  dplyr::filter(status=='singlet') %>%
  ggplot(aes(x=reorder(assignment, -N), y=N)) +
  geom_bar(stat="identity", fill='#FFCCFF')+
  xlab('Genotype assigment')+
  ylab('Cells') +
  labs(subtitle=paste0(gsub('Sample_', '', out),"::", 'Singlets')) +
  geom_text(aes(x=reorder(assignment, -N), y=N, label=N), size=2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  -> Plot_singlets

# Plot 2
clusters_Freq %>% 
  dplyr::filter(status!='singlet') %>%
  ggplot(aes(x=reorder(assignment, -N), y=N)) +
  geom_bar(stat="identity", fill='#CCE5FF')+
  xlab('Genotype assigment')+
  ylab('Cells') +
  labs(subtitle=paste0(gsub('Sample_', '', out),"::", 'Doublets')) +
  geom_text(aes(x=reorder(assignment, -N), y=N, label=N), size=1.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  -> Plot_doublet

print(".. 4. Done plots ..")
#########################################


##########################################
#### Summary Genotyping:  Second part ####
##########################################

#############################
###### Read Cor  table ######
#############################

# Pool directories available
CorTables <- list.files(pattern = "^Pool_Sample*")

# Pool by chr directories available
#CorTables <- CorTables[grepl("chr", CorTables)]

# Pool by chr directories available
CorTables <- CorTables[grepl(Sample, CorTables)]

# Get cor matrix
Cor_by_chr <- Get_GenotypeCor(CorTables)
Cor_by_chr$.id <- gsub(".NA", "", Cor_by_chr$.id)

Cor_by_chr$Cluster <- as.character(Cor_by_chr$Cluster)
#Cor_by_chr$.id <- factor(Cor_by_chr$.id, levels = c(Sample, paste0(Sample,".chr", seq(1,10,1))))

#############################

#########################################################
######       Read bc by chrs and whole genome      ######
#########################################################

Cellgenotypes <- Get_Cellgenotypes(CorTables)

# make list of class by chrs
CellgenotypesDB <- split(Cellgenotypes$assignment, Cellgenotypes$.id)  

# Make  confusion plot cy chr
ConfusionPlot <- lapply(names(CellgenotypesDB)[-c(1)], function(x) Confusiont_plot("P4",x))
names(ConfusionPlot) <- names(CellgenotypesDB)[-c(1)]



#############################
######   Summary Plots  #####
#############################

# Plot 3
ggplot(Cor_by_chr, aes(x = Genotype, y = Cluster, fill = Cor)) +
  geom_tile() +
  scale_fill_viridis(discrete = FALSE) +
  #geom_text(aes(label = Cor), color = "white", size = 1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = paste0("Genotype cor. Sample ",Sample),
       y = "Cell cluster", x = "Genotype") +
  facet_grid(.id ~ .) -> plot_Cor_heatmap

# Plot 4
#Cell_Counts_per_chr$.id <- factor(Cell_Counts_per_chr$.id, levels = paste0("chr", seq(1,10,1)))
ggplot(Cell_Counts_per_chr, aes(y=.id, x=N, fill=assignment)) +
  geom_bar(stat="identity", position=position_dodge())+
  xlab('')+
  ylab('Singlets') +
  labs(subtitle=paste0('Sample ', Sample)) +
  geom_text(aes(label=N), size=1.5, position = position_dodge(width = .9)) +
  scale_fill_viridis(alpha = 0.7, discrete = T)+
  theme_bw() -> Plot_singlets_chrs

pdf(paste0("Plots/SingletSummary_", out, ".pdf"), width = 4, height = 3)
print(Plot_singlets_chrs)
dev.off()

# Plot 5
Plot_1_5 <- {(Plot_singlets/Plot_doublet/Plot_singlets_chrs + 
                plot_layout(heights =  c(0.3, 0.3, 1)))} | plot_Cor_heatmap
# Plot_1_5


pdf(paste0("Plots/Summary_", out, ".pdf"), width = 8, height = 8)
print(Plot_1_5)
dev.off()
print(".. 5. Done saving plots ..")

# Plot 6: confusion matrix
plotlist <- patchwork::wrap_plots(ConfusionPlot, ncol = 3, nrow = 3) 

pdf(paste0("Plots/Confusion_", out, ".pdf"), width = 9, height = 9)
print(plotlist)
dev.off()
