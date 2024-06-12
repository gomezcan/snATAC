
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
#suppressMessages(library(ggpubr))


args <- commandArgs(T)
out <- as.character(args[1])

############################
#### Summary Genotyping ####
############################
# load cluster calls from genotyping
#out <- "Sample_S2_chr10"
clusters <- as.data.table(read.table(paste0('Pool_', out, "/clusters.tsv"), sep = '\t', header = T))


# load cluster assgination
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

#Plot_1 <- {Plot_singlets/Plot_doublet}
Plot_1 <- cowplot::plot_grid(Plot_singlets/Plot_doublet)

pdf(paste0("Plots/SingletSummary_", out, ".pdf"), width = 3, height = 4)
print(Plot_1)
dev.off()

print(".. 5. Done saving plots ..")

write.table(clusters, paste0("Genotyping_Summary_", out, ".txt"), row.names = F, quote = F, sep = '\t')

