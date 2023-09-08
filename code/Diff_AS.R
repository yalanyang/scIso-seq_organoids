
library(tidyr)
library(stringr)
library(dplyr)
library(ggpubr)
library(GenomicRanges)


setwd("/Users/yangyalan/OneDrive - The University of Chicago/Chicago/organiod/Six_organiods_Cellreports/Flair/Neuron_Progenitor")

ES<- read.table("/Users/yangyalan/OneDrive - The University of Chicago/Chicago/organiod/Six_organiods_Cellreports/Flair/diff/Neuron_Progentior.es.fishers.tsv",header=T,sep='\t')
group1 <- "Neuron"
group2 <- "Progenitor"
#ES <- ES %>% filter(ES[,6] <=0.05) 

gene <- read.table('/Users/yangyalan/OneDrive - The University of Chicago/Chicago/organiod/RYA1/Iso-seq/Homo_sapiens.gencode.v40.gene_annotation_table.txt',header=TRUE,sep='\t')
gene1 <- gene %>% select(Gene, GeneSymbol)
ES$Gene <- sapply(strsplit(ES$isoform_ids,split="[_,]"), tail, 1)
ES = ES[,c(2:4,13:14)]
ES =  unname(ES %>% left_join(gene1,by = "Gene"))
ES_merge <- data.frame()
for (i in 1:(nrow(ES)/2)) {
  temp1 <- c(ES[2*i-1,],ES[2*i,])
  ES_merge <-rbind(ES_merge,temp1)
}

ES_merge = ES_merge[,c(1:3,8:12)]
colnames(ES_merge) <- c("coordinate",paste0(group1,"_inclusion"), paste0(group2,"_inclusion"), paste0(group1,"_exclusion"), paste0(group2,"_exclusion"), "Pvalue", "Gene", "Genesymbol")
name1 <- paste0(group1,"_PSI")
name2 <- paste0(group2,"_PSI")
ES_merge <- ES_merge %>% mutate(count= ES_merge[,2]+ES_merge[,3]+ES_merge[,4]+ES_merge[,5], !!name1 := ES_merge[,2]/(ES_merge[,2]+ES_merge[,4]),
                                !!name2 := ES_merge[,3]/(ES_merge[,3]+ES_merge[,5]), 
                                PSI_Diff = ES_merge[,2]/(ES_merge[,2]+ES_merge[,4])- ES_merge[,3]/(ES_merge[,3]+ES_merge[,5]))

ES_merge <- ES_merge %>% separate(coordinate,  c("Chr", "Start","End"), sep = '[-:]')
ES_merge2 <- ES_merge %>% select(Start, End)
##For the exon in s
for (i in 1:nrow(ES_merge)) {
  ES_merge[i,2] <- min (ES_merge2[i,1], ES_merge2[i,2])
  ES_merge[i,3] <- max (ES_merge2[i,1], ES_merge2[i,2])
}
ES_merge$coordinate <- paste0(ES_merge$Chr, ':', ES_merge$Start, "-", ES_merge$End)
ES_merge <- ES_merge %>% filter((Neuron_inclusion+Neuron_exclusion)>10 & (Progenitor_inclusion+Progenitor_exclusion) > 10)


#ES_merge <- ES_merge %>% filter(abs(PSI_Diff)>0.05 & Pvalue <= 0.05) %>% arrange(Pvalue)

write.table(ES_merge, paste0(group1, '_', group2, "_all.alt3.txt"),sep = "\t",quote = F,row.names = F)






