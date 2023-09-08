
library(tidyr)
library(stringr)
library(dplyr)
library(ggpubr)

exon <- read.delim('all.collapse.sorted.filtered.final.multiexon.gtf',header=F)


exon <- exon[grep(pattern="exon",exon$V3),]
exon$V9 <- gsub("transcript_id | gene_id ", "", exon$V9)
exon[, c("Transcript2","Gene")] <- str_split_fixed(exon$V9, ";", 2)
exon$Gene <- gsub(";", "", exon$Gene)
exon$coordination <- paste0(exon$V1, ':', exon$V4, "-", exon$V5,":",exon$V7)
exon = exon %>% select(V1,V4,V5, V7, coordination, Transcript2, Gene)

colnames(exon)[1]<-"Chr"
colnames(exon)[2]<-"Start"
colnames(exon)[3]<-"End"
colnames(exon)[4]<-"Strand"
exon$Transcript = exon$Transcript2

##filter out exons that are overlapped with other exons
#options(bedtools.path = "/Users/yangyalan/miniconda3/bin/")
#library(bedtoolsr)
#overlap_exon <- bt.intersect(a=exon, b=exon,wa=T,wb=T)
#overlap_exon <- overlap_exon%>% dplyr::filter((V1==V9) & ((V2 != V10) | (V3 !=V11)))
#filter_out_exons <- unique(c(overlap_exon$V5, overlap_exon$V13))

#exon <- exon[!(exon$coordination%in%filter_out_exons),] 


##label each exon by gene (exon1...exonN)
exon_label <- unique(exon %>% select(Start, Strand, coordination, Gene))
exon_label$copy <-exon_label$Gene
per_gene_list <- split(x = exon_label[,-5], f =exon_label[,5])
add_exon_grp1 = function(t) {
  if (t$Strand[1] == "+"){
    t <- dplyr::arrange(t,Start)
    t$exon_name = seq(1, nrow(t))
  }
  else{
    t <- dplyr::arrange(t,desc(Start))
    t$exon_name = seq(1, nrow(t))
  }
  return (t)
}
exons_list = lapply(per_gene_list, add_exon_grp1)
exon_label <- do.call(what = rbind, args = exons_list)
#exon_label =  exon_label %>% mutate(exon_name =paste0('exon', exon_name))
exon = exon %>% select(coordination, Transcript2, Gene, Transcript) %>% left_join(exon_label[, c('coordination', 'exon_name')], by="coordination")

##for exon-exon coordination analysis, the first exon and the last exon of each isoform are filtered away. 
per_tran_list <- split(x = exon[,-2], f =exon[,2])
add_exon_grp2 = function(t) {
  t$exon_num = seq(1, nrow(t))
  return (t)
  }
per_gene_list = lapply(per_tran_list, add_exon_grp2)
exon <- do.call(what = rbind, args = per_gene_list)
exon_total = exon %>% group_by(Transcript) %>% dplyr::summarise(n_exon=n())
exon = exon %>% left_join(exon_total, by="Transcript") %>% dplyr::filter(n_exon>2)
exon = exon %>% dplyr::filter(exon_num > 1) %>% filter(exon_num < n_exon) 
length(unique(exon$coordination))

#write.table(exon, "abc.txt",sep = "\t",quote = F,row.names=F)


###count of each exon
count <- read.delim('all.collapse.sorted_RulesFilter_result_classification.Isoform.FL2.txt',header=T)
count <- count[grep(pattern="full-splice_match|novel_not_in_catalog|novel_in_catalog|incomplete-splice_match",count$structural_category),]
count = count %>% select(isoform, associated_gene, FL)
exon_count = exon %>% left_join(count, by=c("Transcript"="isoform"))
exon_count$ID <- paste0(exon_count$coordination,  "-", exon_count$Transcript)


#unique_exon_per_gene <- unique(exon_count %>% select(coordination, exon_name) %>% arrange(exon_name))
exon_count2 = exon_count %>% select(coordination, Gene, Transcript,exon_name)
max_min_exon=exon_count2 %>% group_by(Transcript) %>% dplyr::summarise(max_exon=max(exon_name),min_exon=min(exon_name))
per_gene_list <- split(x = exon_count2[,-2], f =exon_count2[,2])
exon_matrix = function(t) {
  data.2 <-data.frame("coordination"=rep(unique(t$coordination),times=length(unique(t$Transcript))),
                      "Transcript"=rep(unique(t$Transcript),each= length(unique(t$coordination)))
  )
  return (data.2)
}
per_gene = lapply(per_gene_list, exon_matrix)
exon_all <- do.call(what = rbind, args = per_gene)

exon_all$ID <- paste0(exon_all$coordination,  "-", exon_all$Transcript)

exon_anno <- exon_all  %>% left_join(unique(exon_count2[,c("coordination", "exon_name")]), by="coordination") %>% left_join(count[,c("isoform", "associated_gene")], by=c("Transcript"="isoform"))
exon_anno <- exon_anno  %>% left_join(unique(exon_count[,c('Transcript', 'FL')]), by="Transcript")
colnames(exon_anno)[6] <- "out"
exon_anno <- exon_anno  %>% left_join(exon_count[,c('ID', 'FL')], by="ID")
exon_anno$FL[is.na(exon_anno$FL)] <- 0
exon_anno <- transform(exon_anno, out =replace(out, FL>0, "0"))
colnames(exon_anno)[7] <- "inclusion"



## set of alternative exons that met each of the following criteria: 
##(1) ≥10 supporting reads (inclusion + exclusion) in the pseudo-bulk; (2) 0.05 < Ψ < 0.95 at the pseudo-bulk level;
exon_anno[,6:7] <-do.call(cbind,lapply(exon_anno[,6:7],as.numeric))
exon_anno <- exon_anno %>% left_join(max_min_exon, by="Transcript")
exon_anno <- exon_anno%>% dplyr::filter(exon_name >= min_exon & exon_name <= max_exon)
exclusion <- aggregate(exon_anno$out,by=list(type=exon_anno$coordination),sum)
inclusion <- aggregate(exon_anno$inclusion,by=list(type=exon_anno$coordination),sum)
PSI = inclusion  %>% left_join(exclusion, by="type")
colnames(PSI) <- c("coordination","inclusion","exclusion")
PSI <- PSI %>% mutate (psi = inclusion/(inclusion+exclusion), all=inclusion+exclusion)
exon_anno <- exon_anno %>% left_join (PSI[,c('coordination','psi','all')], by="coordination") 
exon_anno <- exon_anno%>% dplyr::filter((psi > 0.05 & psi < 0.95) & all > 20)
exon_anno <- exon_anno[,1:7]


##exon-exon-coordination
exon_total = exon_anno %>% group_by(Transcript) %>% dplyr::summarise(n_exon=n())
exon_anno = exon_anno %>% left_join(exon_total, by="Transcript") %>% dplyr::filter(n_exon>1)

exon_anno$Transcript2 = exon_anno$Transcript
per_tran_list <- split(x = exon_anno[,-9], f =exon_anno[,9])

coexon = function(t) {
  tmp<-matrix(0,nrow=(nrow(t)*(nrow(t)-1)/2),ncol=9)
  k=0        
  for (i in 1:(nrow(t)-1)){
    for (j in 2:nrow(t)){
    if (i<j){
      k <-k+1;
    tmp[k,1] = paste0(t$coordination[i],'-', t$coordination[j]);
    tmp[k,2] = t$associated_gene[1];
    tmp[k,3] = t$Transcript[1];
    tmp[k,4] = t$inclusion[i];
    tmp[k,5] = t$inclusion[j];
    tmp[k,8] = t$coordination[i]
    tmp[k,9] = t$coordination[j]
  if(t$inclusion[i] > 0 & t$inclusion[j] > 0){
   tmp[k,6] = t$inclusion[i]; 
  tmp[k,7] <- "in-in";
  }
  else if(t$inclusion[i] > 0 & t$inclusion[j] == 0){
      tmp[k,6] = t$inclusion[i];
      tmp[k,7] <- "in-out";
     }
      else if (tmp[k,4] == 0 & tmp[k,5] > 0){
      tmp[k,6] = t$inclusion[j];
      tmp[k,7] <- "out-in";
      }else if (tmp[k,4] == 0 & tmp[k,5] == 0){
      tmp[k,6] = t$out[i];
      tmp[k,7] <- "out-out";
      }
      }
    }
  }
  return(tmp)
}
coexon_list = lapply(per_tran_list, coexon)
co_exon <- data.frame(do.call(what = rbind, args = coexon_list))
colnames(co_exon) <- c("exonID","Gene","Transcript","1stexon_count","2ndexon_count","count","type","exon1","exon2")



#transfer to matrix
co_exon$exonID2 = co_exon$exonID
per_coexon_list <- split(x = co_exon[,-10], f =co_exon[,10])
type_count = function(t) {
  tmp<-matrix(0,nrow=1,ncol=8)
  tmp[1,1] <- t$exonID[1];
  tmp[1,2] <- t$Gene[1];
  tmp[1,7] = t$exon1[1]
  tmp[1,8] = t$exon2[1]
  tmp[1,3] <- sum(as.numeric(t[grep(pattern = "in-in",t$type),]$count))
  tmp[1,4] <- sum(as.numeric(t[grep(pattern = "in-out",t$type),]$count)) 
 tmp[1,5] <- sum(as.numeric(t[grep(pattern = "out-in",t$type),]$count)) 
  tmp[1,6] <- sum(as.numeric(t[grep(pattern = "out-out",t$type),]$count)) 
  return (tmp)
}
coexon_list = lapply(per_coexon_list, type_count)
co_exon2 <- as.data.frame(do.call(what = rbind, args = coexon_list))

co_exon2[,3:6] <-do.call(cbind,lapply(co_exon2[,3:6],as.numeric))
colnames(co_exon2) <- c("exonID","Gene","in_in","in_out","out_in","out_out","exon1","exon2")
exon_order_gene <- unique(exon[c("coordination","exon_name")])
co_exon2 <- co_exon2  %>% left_join(exon_order_gene, by=c("exon1"="coordination")) %>% left_join(exon_order_gene, by=c("exon2"="coordination"))
co_exon2 <- co_exon2 %>% mutate(sum=in_in+in_out+out_in+out_out,OD=log2((in_in+1)*(out_out+1)/((in_out+1)*(out_in+1))))

co_exon2$exon1_copy <- co_exon2$exon1
co_exon2$exon2_copy <- co_exon2$exon2
co_exon2 <- co_exon2 %>% separate(exon1_copy,  c(NA, "exon1_start", "exon1_end"), sep = '[:-]')%>% separate(exon2_copy,  c(NA, "exon2_start", "exon2_end"), sep = '[:-]')


co_exon_filter_away <- co_exon2 %>% dplyr::filter((exon2_start>=exon1_start & exon2_start <= exon1_end)|(exon2_end >=exon1_start & exon2_end <= exon1_end))
co_exon_filter_away <- as.matrix(sort(unique(c(co_exon_filter_away$exon1,co_exon_filter_away$exon2))))
co_exon2 <- co_exon2[!(co_exon2$exon1%in%co_exon_filter_away[,1]),] 
co_exon2 <- co_exon2[!(co_exon2$exon2%in%co_exon_filter_away[,1]),] 
co_exon2 <- co_exon2 %>% dplyr::filter(exon2_start > exon1_end | exon2_end < exon1_start)


chisq_test = function(exonID){
  Xsq <- chisq.test(matrix(as.numeric(unlist(c(co_exon2[co_exon2$exonID==exonID, ][3:6]))),nrow=2,ncol=2), correct=F)
  return(Xsq$p.value)
}
coexon_list2 = lapply(co_exon2$exonID, chisq_test)
co_exon3 <- data.frame(do.call(what = rbind, args = coexon_list2))
co_exon3 <-cbind(co_exon2,co_exon3)
colnames(co_exon3)[17]<-"Pvalue"
co_exon3 <- co_exon3 %>% dplyr::filter(Pvalue != "NaN")
co_exon3$FDR <-  p.adjust(co_exon3$Pvalue, method = "BY", n = length(co_exon3$Pvalue))

co_exon3 <- co_exon3 %>% dplyr::filter(FDR < 0.001 & abs(OD)>2) %>% arrange(FDR)
write.table(co_exon3, "co_exon.txt",sep = "\t",quote = F,row.names=F)

