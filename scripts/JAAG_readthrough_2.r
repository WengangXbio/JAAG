.libPaths('~/schoenebeck_group/WENGANG/R_lib/')
library(splitstackshape)
library(dplyr)
library(TCseq)
library(purrr)
library(cluster)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
reassign=read.table(args[1],head=F)
bed=read.table(args[2],head=F)
blast=read.table(args[3],head=F)
rm_list=read.table(args[4],head=F)
rm_list=data.frame(rm_list,ind=1)
colnames(rm_list)=c("V2","ind")
blast =left_join(blast,rm_list,by="V2")
blast =blast[-which(blast[,13]==1),]
blast=blast[,c(1:12)]
blast[,1] <- cSplit(blast, 'V1', sep="(", type.convert=FALSE)[,12]
blast[,13] <- cSplit(blast, 'V1', sep=";", type.convert=FALSE)[,12]
coding_ls=data.frame(V4=unique(blast[,1]),ind=1)
coding_bed=left_join(bed,coding_ls,by="V4")
coding_bed=coding_bed[which(coding_bed[,13]==1),c(1:12)]
rethgene <- data.frame(unique(cSplit(reassign, 'V1', sep=";", type.convert=FALSE)[,2]),ind=1)
colnames(rethgene) = c('V4_1','ind')
coding_bed[,13] <-cSplit(coding_bed, 'V4', sep=";", type.convert=FALSE)[,12]
coding_bed=left_join(coding_bed,rethgene,by="V4_1")
coding_bed_readthrougth=coding_bed[!is.na(coding_bed[,14]),c(1:12)]
coding_bed_normal=coding_bed[is.na(coding_bed[,14]),c(1:12)]
colnames(reassign) = c('V4','new')
coding_bed_readthrougth=left_join(coding_bed_readthrougth,reassign,by="V4")
coding_bed_readthrougth=na.omit(coding_bed_readthrougth)
coding_bed_readthrougth[,4]=coding_bed_readthrougth[,13]
coding_bed_readthrougth=coding_bed_readthrougth[,-13]
output_bed=rbind(coding_bed_normal,coding_bed_readthrougth)
write.table(output_bed,'merged_rmpolya_rmsingle_coding.readthrough_reassign.bed',quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
