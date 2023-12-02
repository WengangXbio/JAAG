forge_pm <- function(asda){
rem=matrix(0,nrow=length(unique(asda[,1])),ncol=length(unique(asda[,2])))
colnames(rem)=unique(asda[,2])
rownames(rem)=unique(asda[,1])
for(nd in 1:nrow(asda)){
rem[which(rownames(rem)==asda[nd,1]),which(colnames(rem)==asda[nd,2])]=1
}
rem
}
Average_Silhouette  <- function (ttput) {
zhi=qr(ttput)$rank
avg_sil <- function(k) {
  km.res <- kmeans(ttput, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(ttput))
  mean(ss[, 3])
}
# Compute and plot wss for k = 2 to k = 15
k.values <- 2:(zhi-1)
if(zhi >5)
k.values <- 2:5
# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)
which(avg_sil_values==max(avg_sil_values))+1
}
forge_group <- function(aq,n1){
test_group1=list()
p=1
for(i in 1:n1){
if(length(names(aq$cluster[which(aq$cluster==i)]))==1)
next
test_group1[[p]]=names(aq$cluster[which(aq$cluster==i)])
p=p+1
}
test_group1
}
cont_super <- function(t1,t2,group,gene){
g1=data.frame(V4=group[[t1]],ind=1)
v1=left_join(gene,g1,by="V4")
g2=data.frame(V4=group[[t2]],ind=1)
v2=left_join(gene,g2,by="V4")
v1=na.omit(v1)
v2=na.omit(v2)
v1m=data.table(chr=v1[,1],start=as.numeric(v1[,2]),end=as.numeric(v1[,3]),id=v1[,4])
v2m=data.table(chr=v2[,1],start=as.numeric(v2[,2]),end=as.numeric(v2[,3]),id=v2[,4])
setkey(v1m,chr, start, end)
setkey(v2m,chr, start, end)
overlap=foverlaps(v1m,v2m, nomatch = 0)
if(nrow(overlap)!=0){
pair=data.frame(overlap$id,overlap$i.id)
rid=distinct(pair)
tid=data.frame(table(rid[,2]))
tid=tid[order(-tid$Freq),]
rm_trans_n=ceiling(length(group[[t1]])*1/3)
rm_trans=as.character(tid[,1][1:rm_trans_n])
rm_trans=data.frame(i.id=rm_trans,ind=1)
overlap=left_join(overlap,rm_trans,by="i.id")
overlap=overlap[-which(overlap$ind==1),]
}
length(unique(overlap$id))/length(group[[t2]])
}

judge_fixed_group <- function(judge_mat_1,test_group_1){
p_order=c()
for(hg in 1:nrow(judge_mat_1)){
p_order=c(p_order,sum(judge_mat_1[,hg],judge_mat_1[hg,]))
}
p_order=data.frame(cluster=c(1:nrow(judge_mat_1)),score=p_order)
p_order=p_order[order(-p_order$score),]
for(ds in 1:nrow(p_order)){
p_order[ds,2]=10^ds
}
a=c(1:length(test_group_1))
poss_com=list()
for(wsd in 2:length(a)){
poss_com=append (poss_com,combn(a,wsd,simplify=F))
}
hs=c()
for(sdw in 1:length(poss_com)){
cr=unlist(poss_com[sdw])
if (max(judge_mat_1[cr,cr])<=0.34)
hs=c(hs,sdw)
}
score=c()
for(is in 1:length(hs) ){
mm=data.frame(cluster=unlist(poss_com[[hs[is]]]),1)
score=c(score,sum(na.omit(left_join(p_order,mm,by='cluster'))$score))
}
test_group_1[poss_com[[hs[which(score==max(score))]]]]
}
remove_hybird <- function(test_gene_2,fix_group_2){
return_group=list()
for(xs in 1:length(fix_group_2)){
dsa=unlist(fix_group_2[xs])
grouped_1=data.frame(V4=dsa,ind=1)
cross=left_join(test_gene_2,grouped_1,by='V4')
group_bed=na.omit(cross)
other_bed=cross[-which(cross$ind==1),]
subject_1=data.table(chr=group_bed[,1],start=as.numeric(group_bed[,2]),end=as.numeric(group_bed[,3]),id=group_bed[,4])
query_1=data.table(chr=other_bed[,1],start=as.numeric(other_bed[,2]),end=as.numeric(other_bed[,3]),id=other_bed[,4])
setkey(subject_1,chr, start, end)
setkey(query_1,chr, start, end)
overlap=foverlaps(query_1,subject_1, nomatch = 0)
if(nrow(overlap)!=0){
pair=data.frame(overlap$id,overlap$i.id)
rid=distinct(pair)
tid=data.frame(table(rid[,1]))
tid=tid[order(-tid$Freq),]
rm_trans_n=ceiling(length(dsa)*1/3)
rm_trans=as.character(tid[,1])[c(1:rm_trans_n)]
keeptrans=setdiff(dsa,rm_trans)
return_group[[xs]]=keeptrans
}else{
return_group[[xs]]=dsa
}
}
return_group
}
judge_transcripts <-function(test_gene_1,fix_group_1){
trans_name=unique(test_gene_1[,4])
final_group=list()
final_group[c(1:(length(fix_group_1)+2))]='seed'
for(sq in 1:length(trans_name)){
tested=test_gene_1[which(test_gene_1[,4]==trans_name[sq]),]
query=data.table(chr=tested[,1],start=as.numeric(tested[,2]),end=as.numeric(tested[,3]))
setkey(query,chr, start, end)
cvz=c()
for(gg in 1:length(fix_group_1)){
grouped=data.frame(V4=unlist(fix_group_1[gg]),ind=1)
fix_group_bed=na.omit(left_join(test_gene_1,grouped,by='V4'))
subject=data.table(chr=fix_group_bed[,1],start=as.numeric(fix_group_bed[,2]),end=as.numeric(fix_group_bed[,3]))
setkey(subject,chr, start, end)
overlap=foverlaps(query,subject, nomatch = 0)
if(nrow(overlap))
cvz=c(cvz,gg)
}
if(length(cvz)==1){
final_group[[cvz]]=c(final_group[[cvz]],as.character(trans_name[sq]))
}else if (length(cvz)>1){
final_group[[(length(fix_group_1)+1)]]=c(final_group[[(length(fix_group_1)+1)]],as.character(trans_name[sq]))
}else if (length(cvz)==0){
final_group[[(length(fix_group_1)+2)]]=c(final_group[[(length(fix_group_1)+2)]],as.character(trans_name[sq]))
}
}
lapply(final_group, function(x) x[-1])
}

stat <- function(output2){
statreturn=data.frame(V1=0,V2=0,V3=0,V4=0,V5=0)
for(u in 1:length(output2)){
toutput2=output2[[u]]
na=toutput2[[length(toutput2)]]
reth=toutput2[[length(toutput2)-1]]
lg=c()
for(hd in 1:(length(toutput2)-2)){
lg=c(lg,length(toutput2[[hd]]))
}
addstat=data.frame(V1=length(na),V2=length(reth),V3=min(lg),V4=sum(lg),V5=length(lg))
statreturn=rbind(statreturn,addstat)
}
colnames(statreturn)=c('NA','readthrough','min_count','sum_count','n_group')
statreturn[-1,]
}

clean_up <- function(output3){
ggg=data.frame(V1=0,V2=0)
for(u in 1:length(output3)){
toutput3=output3[[u]]
lg=c()
for(hd in 1:(length(toutput3)-2)){
if(length(toutput3[[hd]])<=1)
next
name=data.frame(V1=toutput3[[hd]])
name_new=cSplit(name, 'V1', sep=";", type.convert=FALSE)
name_new=cSplit(name_new, 'V1_2', sep=".", type.convert=FALSE)
name_new$dash='_'
name_new$dot='.'
name_new$semi=';'
name_new$hd=hd
cols=c('V1_1','dash','hd','semi','V1_1','dash','hd','dot','V1_2_2')
name$V2 <- apply(name_new[,c(1,4,7,6,1,4,7,5,3)],1,paste,collapse = "" )
ggg=rbind(ggg,name)
}
}
ggg[-1,]
}

clean_up_readthrough <- function(output3){
ggg=data.frame(V1=0,V2=0)
for(u in 1:length(output3)){
toutput3=output3[[u]]
gh=length(toutput3)
name=data.frame(V1=toutput3[[gh-1]])
if(nrow(name)==0)
next
name_new=cSplit(name, 'V1', sep=";", type.convert=FALSE)
name_new=cSplit(name_new, 'V1_2', sep=".", type.convert=FALSE)
name_new$dash='_'
name_new$dot='.'
name_new$semi=';'
name_new$hd='readthrough'
cols=c('V1_1','dash','hd','semi','V1_1','dash','hd','dot','V1_2_2')
name$V2 <- apply(name_new[,c(1,4,7,6,1,4,7,5,3)],1,paste,collapse = "" )
ggg=rbind(ggg,name)
}
ggg[-1,]
}


delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
    x.list[unlist(lapply(x.list, length) != 0)]
}

.libPaths('~/schoenebeck_group/WENGANG/R_lib/')
library(splitstackshape)
library(dplyr)
library(TCseq)
library(purrr)
library(cluster)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
genebed=read.table(args[1],head=F)
blast=read.table(args[2],head=F)
rm_list=read.table(args[3],head=F)
rm_list=data.frame(rm_list,ind=1)
colnames(rm_list)=c("V2","ind")
blast =left_join(blast,rm_list,by="V2")
blast =blast[-which(blast[,13]==1),]
blast=blast[,c(1:12)]
blast[,1] <- cSplit(blast, 'V1', sep="(", type.convert=FALSE)[,12]
blast[,13] <- cSplit(blast, 'V1', sep=";", type.convert=FALSE)[,12]
genebed[,7] <- cSplit(genebed, 'V4', sep=";", type.convert=FALSE)[,6]
genen=unique(blast[,13])
output=list()
g=1
for(f in 1:length(genen)){
subblast=blast[which(blast[,13]==genen[f]),]
pm <- forge_pm(subblast) 
if(qr(pm)$rank <2||dim(pm)[1]<5)
next
if(qr(pm)$rank==2){
n <- 2
}else{
n <- Average_Silhouette(pm)[1]
}
km.res <- kmeans(pm, centers = n, nstart = 25)
test_group <- forge_group(km.res,n)
if(length(test_group)==1)
next
judge_mat <- matrix(0,ncol=length(test_group),nrow=length(test_group))
for(j in 1:(length(test_group)-1)){
for(k in (j+1):length(test_group)){
judge_mat[j,k]=min(cont_super(j,k,test_group,genebed),cont_super(k,j,test_group,genebed))
judge_mat[k,j]=min(cont_super(j,k,test_group,genebed),cont_super(k,j,test_group,genebed))
}
}
diag(judge_mat) <- NA
if(min(judge_mat, na.rm = TRUE)>0.34)
next
diag(judge_mat) <- 0
fix_group <- judge_fixed_group(judge_mat,test_group)
if(length(fix_group)<2)
next
pure_fix_group <- remove_hybird(genebed,fix_group)
lsoutput <- judge_transcripts(genebed[which(genebed[,7]==genen[f]),],pure_fix_group)
osubblast=data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0)
colnames(osubblast)=colnames(subblast)
for(o in 1:length(unique(subblast[,1]))){
ssblast=subblast[which(subblast[,1]==unique(subblast[,1])[o]),]
addssblast=ssblast[which(ssblast[,12]>=max(ssblast[,12])*0.95),]
osubblast=rbind(osubblast,addssblast)
}
osubblast=osubblast[-1,]
for(kk in 1:(length(lsoutput)-3)){
for(zz in (kk+1):(length(lsoutput)-2)){
if(length(lsoutput[[zz]])*length(lsoutput[[kk]])==0)
next
gg1=data.frame(V1=lsoutput[[kk]],ind=1)
gg2=data.frame(V1=lsoutput[[zz]],ind=1)
c1=na.omit(left_join(osubblast,gg1,by="V1"))
c2=na.omit(left_join(osubblast,gg2,by="V1"))
if(length(intersect(unique(c1[,2]),unique(c2[,2])))!=0){
lsoutput[[kk]]=c(lsoutput[[zz]],lsoutput[[kk]])
lsoutput[[zz]]=character(0)
}
}
}
youtput <- delete.NULLs(lsoutput[1:(length(lsoutput)-2)])
youtput[length(youtput)+1] <- lsoutput[(length(lsoutput)-1)]
youtput[length(youtput)+1] <- lsoutput[length(lsoutput)]
if(length(youtput)<=3)
next
youtput2 <- judge_transcripts(genebed[which(genebed[,7]==genen[f]),],youtput[c(1:(length(youtput)-2))])
output[[g]] <- youtput2
g=g+1
}
readthrough_stat=stat(output)
output_rm=output[- which(readthrough_stat[,2]>=readthrough_stat[,4])]
rename <- clean_up(output_rm)
rename_readthrough <- clean_up_readthrough(output_rm)

###second round###
output_rmout=output[which(readthrough_stat[,2]>=readthrough_stat[,4])]
sename=as.data.frame(unlist(output_rmout))
colnames(sename)='V1'
sename <- cSplit(sename, 'V1', sep=";", type.convert=FALSE)[,1]
sename =unlist(c(unique(sename)))
gg=1
output3=list()
for(l in 1:length(sename)){
subblast=blast[which(blast[,13]==sename[l]),]
pm <- forge_pm(subblast) 
if(qr(pm)$rank <3||dim(pm)[1]<5)
next
n <- Average_Silhouette(pm)[1]
km.res <- kmeans(pm, centers = n, nstart = 25)
test_group <- forge_group(km.res,n)
if(length(test_group)==1)
next
judge_mat <- matrix(0,ncol=length(test_group),nrow=length(test_group))
for(j in 1:(length(test_group)-1)){
for(k in (j+1):length(test_group)){
judge_mat[j,k]=min(cont_super(j,k,test_group,genebed),cont_super(k,j,test_group,genebed))
judge_mat[k,j]=min(cont_super(j,k,test_group,genebed),cont_super(k,j,test_group,genebed))
}
}
diag(judge_mat) <- NA
if(min(judge_mat, na.rm = TRUE)>0.34)
next
diag(judge_mat) <- 0
fix_group <- judge_fixed_group(judge_mat,test_group)
if(length(fix_group)<2)
next
pure_fix_group <- remove_hybird(genebed,fix_group)
lsoutput <- judge_transcripts(genebed[which(genebed[,7]==sename[l]),],pure_fix_group)
osubblast=data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0)
colnames(osubblast)=colnames(subblast)
for(o in 1:length(unique(subblast[,1]))){
ssblast=subblast[which(subblast[,1]==unique(subblast[,1])[o]),]
addssblast=ssblast[which(ssblast[,12]>=max(ssblast[,12])*0.95),]
osubblast=rbind(osubblast,addssblast)
}
osubblast=osubblast[-1,]
for(kk in 1:(length(lsoutput)-3)){
for(zz in (kk+1):(length(lsoutput)-2)){
if(length(lsoutput[[zz]])*length(lsoutput[[kk]])==0)
next
gg1=data.frame(V1=lsoutput[[kk]],ind=1)
gg2=data.frame(V1=lsoutput[[zz]],ind=1)
c1=na.omit(left_join(osubblast,gg1,by="V1"))
c2=na.omit(left_join(osubblast,gg2,by="V1"))
if(length(intersect(unique(c1[,2]),unique(c2[,2])))!=0){
lsoutput[[kk]]=c(lsoutput[[zz]],lsoutput[[kk]])
lsoutput[[zz]]=character(0)
}
}
}
youtput <- delete.NULLs(lsoutput[1:(length(lsoutput)-2)])
youtput[length(youtput)+1] <- lsoutput[(length(lsoutput)-1)]
youtput[length(youtput)+1] <- lsoutput[length(lsoutput)]
if(length(youtput)<=3)
next
youtput2 <- judge_transcripts(genebed[which(genebed[,7]==sename[l]),],youtput[c(1:(length(youtput)-2))])
youtput3 <- judge_transcripts(genebed[which(genebed[,7]==sename[l]),],youtput2[c(1:(length(youtput2)-2))])
output3[[gg]] <- youtput2
gg=gg+1
}
readthrough_stat3=stat(output3)
output3_rm=output3[- which(readthrough_stat3[,2]>=readthrough_stat3[,4])]
rename3 <- clean_up(output3_rm)
rename_readthrough3 <- clean_up_readthrough(output3_rm)
if(length(rename3)!=0){
renames <- rbind(rename,rename3)
rename_readthroughs <- rbind(rename_readthrough,rename_readthrough3)
}
write.table(renames,'merged_rmpolya_rmsingle_coding.readthrough.reassign',quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(rename_readthroughs,'merged_rmpolya_rmsingle_coding.readthrough.collection',quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
