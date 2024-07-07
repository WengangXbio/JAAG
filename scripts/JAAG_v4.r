merge_fun_p <- function(inputbed,inputibed,blast,geneid,jet_f,distance_f){
for(i in 1:nrow(inputbed)){  #for(i in 333:350){
subblast=blast[which(blast[,1]==as.character(inputbed[i,4])),]
if(nrow(subblast)==0)
next
j=c((i+1):(i+4))
query=data.frame(sub=as.character(inputbed[j,4]),index=1)
mda=left_join(blast,query,by="sub")
mda=na.omit(mda)
subg=as.character(unique(intersect(mda$gene,subblast$gene)))
if(!length(subg))
next
yeh=data.frame(gene=subg,index=1)
candilist=unique(na.omit(left_join(mda,yeh,by="gene"))[,1])
dis=c(inputibed[which(inputibed[,4]==candilist[1]),2]-inputibed[which(inputibed[,4]==as.character(inputbed[i,4])),3],inputibed[which(inputibed[,4]==candilist[2]),2]-inputibed[which(inputibed[,4]==as.character(inputbed[i,4])),3])
if(max(dis)>distance_f)
candilist=candilist[-which(dis>distance_f)]
if(length(candilist)){
candilist2=c()
for(u in 1:length(candilist)){
submda=mda[which(mda$sub==candilist[u]),]
if(length(intersect(submda$gene,subblast$gene)))
candilist2=c(candilist2,candilist[u])
}
if(length(candilist2)){
for(d in 1:length(candilist2)){
submda=mda[which(mda$sub==candilist2[d]),]
judgegene=intersect(submda$gene,subblast$gene)[1]
cv1=subblast[which(subblast$gene==judgegene),]
cv2=submda[which(submda$gene==judgegene),]
if((max(cv1$send)-jet_f)<=min(cv2$sgene)){
if(geneid[which(geneid[,1]==subblast$sub[1]),2]=='NA'){
if(geneid[which(geneid[,1]==submda$sub[1]),2]=='NA'){
geneid[which(geneid[,1]==subblast$sub[1]),2]=judgegene[1] #paste("MERGE",dd,sep='')  #subblast$gene[1]
}else{
geneid[which(geneid[,1]==subblast$sub[1]),2]=geneid[which(geneid[,1]==submda$sub[1]),2]
}
}
geneid[which(geneid[,1]==submda$sub[1]),2]=geneid[which(geneid[,1]==subblast$sub[1]),2]
#dd=dd+1
}
}
}
}
}
geneid
}

merge_fun_n <- function(inputbed,inputibed,blast,geneid,jet_f,distance_f){
for(i in 1:nrow(inputbed)){  
subblast=blast[which(blast[,1]==as.character(inputbed[i,4])),]
if(nrow(subblast)==0)
next
j=c((i+1):(i+4))
query=data.frame(sub=as.character(inputbed[j,4]),index=1)
mda=left_join(blast,query,by="sub")
mda=na.omit(mda)
subg=as.character(unique(intersect(mda$gene,subblast$gene)))
if(!length(subg))
next
yeh=data.frame(gene=subg,index=1)
candilist=unique(na.omit(left_join(mda,yeh,by="gene"))[,1])
dis=c(inputibed[which(inputibed[,4]==candilist[1]),2]-inputibed[which(inputibed[,4]==as.character(inputbed[i,4])),3],inputibed[which(inputibed[,4]==candilist[2]),2]-inputibed[which(inputibed[,4]==as.character(inputbed[i,4])),3])
if(max(dis)>distance_f)
candilist=candilist[-which(dis>distance_f)]
if(length(candilist)){
candilist2=c()
for(u in 1:length(candilist)){
submda=mda[which(mda$sub==candilist[u]),]
if(length(intersect(submda$gene,subblast$gene)))
candilist2=c(candilist2,candilist[u])
}
if(length(candilist2)){
for(d in 1:length(candilist2)){
submda=mda[which(mda$sub==candilist2[d]),]
judgegene=intersect(submda$gene,subblast$gene)[1]
cv1=subblast[which(subblast$gene==judgegene),]
cv2=submda[which(submda$gene==judgegene),]
if((min(cv1$sgene)+jet_f)>=max(cv2$send)){
if(geneid[which(geneid[,1]==subblast$sub[1]),2]=='NA'){
if(geneid[which(geneid[,1]==submda$sub[1]),2]=='NA'){
geneid[which(geneid[,1]==subblast$sub[1]),2]=judgegene[1] #paste("MERGE",dd,sep='')  #subblast$gene[1]
}else{
geneid[which(geneid[,1]==subblast$sub[1]),2]=geneid[which(geneid[,1]==submda$sub[1]),2]
}
}
geneid[which(geneid[,1]==submda$sub[1]),2]=geneid[which(geneid[,1]==subblast$sub[1]),2]
#dd=dd+1
}
}
}
}
}
geneid
}

bestHSPs <- function(dogblast1,SIM_f){
dogblastname=unique(as.matrix(dogblast1)[,1])
newdogblast1=data.frame(0,0,0,0,0,0,0)
colnames(newdogblast1)=c("sub","gene","qstart","qend","sgene","send","score")
for(k in 1:length(dogblastname)){ 
subblast=dogblast1[which(dogblast1[,1]==as.character(dogblastname[k])),]
v1 = distinct(subblast,gene,.keep_all = TRUE)
msubblast = anti_join(subblast,v1,by = c("sub", "gene", "qstart", "qend", "sgene", "send", "score"))
if(nrow(msubblast)!=0){
for(m in 1:nrow(msubblast)){
v2=msubblast[m,]
v1ma=data.table(chr=v1$gene,start=as.numeric(v1$qstart),end=as.numeric(v1$qend))
v2ma=data.table(chr=v2$gene,start=as.numeric(v2$qstart),end=as.numeric(v2$qend))
setkey(v1ma,chr, start, end)
setkey(v2ma,chr, start, end)
overlap1=foverlaps(v1ma,v2ma, nomatch = 0)
v1mb=data.table(chr=v1$gene,start=as.numeric(v1$sgene),end=as.numeric(v1$send))
v2mb=data.table(chr=v2$gene,start=as.numeric(v2$sgene),end=as.numeric(v2$send))
setkey(v1mb,chr, start, end)
setkey(v2mb,chr, start, end)
overlap2=foverlaps(v1mb,v2mb, nomatch = 0)
if((nrow(overlap1)+nrow(overlap2))!=0){
mover2 <- data.table(
    mgene = overlap1$geneid,
    start = overlap1[, ifelse(start > i.start, start, i.start)],
    end = overlap1[, ifelse(end < i.end, end, i.end)])
identi1=sum(mover2$end-mover2$start)/3  #nucleitide overlap
mover3 <- data.table(
    mgene = overlap2$geneid,
    start = overlap2[, ifelse(start > i.start, start, i.start)],
    end = overlap2[, ifelse(end < i.end, end, i.end)])
identi2=sum(mover3$end-mover3$start)    #peptide overlap
if((identi1+identi2)<30)
v1=rbind(v1,v2)
}else{
v1=rbind(v1,v2)
}
}
}
vv1=v1[,c(2,7)]
besth=aggregate(. ~ gene, data=vv1, sum)
vc=besth[which(besth$score>=max(besth$score)*SIM_f),1]
edi=data.frame(gene=vc,index=1)
addblast=left_join(v1,edi,by="gene")
addblast=na.omit(addblast)[,-8]
newdogblast1=rbind(newdogblast1,addblast)
}
newdogblast1[-1,]
}

last_third_clean <- function(candi_bed61,merge_matrix1,distance_f){
nm=unique(merge_matrix1[,2])
new_merge_matrix1=data.frame(0,0)
colnames(new_merge_matrix1)=c("V4","id")
for(r in 1:length(nm)){
fy=unique(candi_bed61[which(candi_bed61[,8]==as.character(nm[r])),4])
fu=candi_bed61[which(candi_bed61[,8]==as.character(nm[r])),]
chr=c()
star=c()
for(v in 1:length(fy)){
subfu=fu[which(fu[,4]==fy[v]),]
chr=c(chr,subfu[1,1])
star=c(star,min(subfu[,2]))
}
fff=data.frame(chr,star,fy)
fff=fff[with(fff, order(chr, star)), ]
hzlabel=1
label=1
if(length(fy)==1)
next
for(s in 2:length(fy)){
hzchr=fff[s-1,1]
hzstar=fff[s-1,2]
if(hzchr==fff[s,1]&abs(hzstar-fff[s,2])<=distance_f){
hzlabel=c(hzlabel,label)
}else if(hzchr==fff[s,1]&abs(hzstar-fff[s,2])>distance_f){
label=label+1
hzlabel=c(hzlabel,label)
}else if(hzchr!=fff[s,1]){
label=label+1
hzlabel=c(hzlabel,label)
}
}
fff=data.frame(fff,hzlabel,id=paste("merge.",r,sep=""))
addfff=data.frame(V4=fff[,3],id=as.matrix(apply(fff[,c(5,4)],1,paste,collapse=";")))
new_merge_matrix1=rbind(new_merge_matrix1,addfff)
}
new_merge_matrix1[-1,]
}
last_second_clean <- function(new_merge_matrix1){
new2_merge_matrix =cSplit(new_merge_matrix1, 'V4', sep="_", type.convert=FALSE)
new2_merge_matrix=na.omit(new2_merge_matrix)[,1:2]
rmfu=unlist(new2_merge_matrix[duplicated(new2_merge_matrix),][,1])
rmfu= data.frame(id=rmfu,ind=1)
new_merge_matrix1=left_join(new_merge_matrix1,rmfu,by='id')
new_merge_matrix1[-which(new_merge_matrix1[,3]==1),-3]
}

library(splitstackshape)
library(dplyr)
library(data.table)
library(TCseq)
options("scipen"=15)
SIM=0.75 #similarity
jet=30   #error tolerant
distance=500000
args <- commandArgs(trailingOnly = TRUE)
ibed=read.table(args[1],head=F,sep='\t')
dogblast=read.table(args[3],head=F,sep='\t')
rm_list=read.table(args[4],head=F)
dogblast=dogblast[,c(1,2,7,8,9,10,12)]
colnames(dogblast)=c("sub","gene","qstart","qend","sgene","send","score")
rm_list=data.frame(rm_list,ind=1)
colnames(rm_list)=c("gene","ind")
dogblast =left_join(dogblast,rm_list,by="gene")
dogblast =dogblast[-which(dogblast[,8]==1),]
dogblast=dogblast[,c(1:7)]

matchid=unique(dogblast[,1])
setmatrix=data.frame(matchid,ind=1)
setmatrix=cSplit(setmatrix, 'matchid', sep="(", type.convert=FALSE)
setmatrix=data.frame(V4=setmatrix[,2],ind=1)
colnames(setmatrix)=c("V4","ind")
setmatrix=distinct(setmatrix)
ibed=left_join(ibed,setmatrix,by="V4")
ibed=na.omit(ibed)
ibedp=ibed[which(ibed[,6]=='+'),]
ibedn=ibed[which(ibed[,6]=='-'),]
ibedp=ibedp[with(ibedp, order(ibedp[,1], ibedp[,2])), ]
ibedn=ibedn[with(ibedn, order(ibedn[,1], ibedn[,3])), ]
dogblast=cSplit(dogblast, 'sub', sep="(", type.convert=FALSE)
dogblast=dogblast[,c(7,1:6)]
colnames(dogblast)=c("sub","gene","qstart","qend","sgene","send","score")
newdogblast = bestHSPs(dogblast,SIM)

geneid6=as.matrix(data.frame(name=ibed[,4],id='NA'))
geneid6=merge_fun_p(ibedp,ibed,newdogblast,geneid6,jet,distance)
geneid6=merge_fun_n(ibedn,ibed,newdogblast,geneid6,jet,distance)


geneid6=geneid6[,c(1,2)]
merge_matrix=geneid6[which(geneid6[,2]!='NA'),]
candim=data.frame(V4=geneid6[,1],ind=1)
colnames(merge_matrix)=c("V4","V7")
merge_matrix=as.data.frame(merge_matrix)

bed6=read.table(args[2],head=F,sep='\t')
candi_bed6=left_join(bed6,candim,by="V4")
candi_bed6=na.omit(candi_bed6)
candi_bed6=left_join(candi_bed6,merge_matrix,by="V4")

new_merge_matrix = last_third_clean(candi_bed6,merge_matrix,distance)
#new_merge_matrix = last_second_clean(new_merge_matrix)

candi_bed6=left_join(candi_bed6,new_merge_matrix,by="V4")
nm=unique(new_merge_matrix[,2])
querybed=c()
subjectbed=c()
que=1
for(r in 1:length(nm)){
fy=unique(candi_bed6[which(candi_bed6[,9]==as.character(nm[r])),4])
querybed=c(querybed,fy)
mg=paste("Merge.",que,"-T",sep="")
candi_bed6[which(candi_bed6[,9]==nm[r]),4]=mg
subjectbed=c(subjectbed,rep(mg,length(fy)))
que=que+1
}
clue=data.frame(querybed,subjectbed)
candi_bed6=candi_bed6[,c(1:6)]

um=unique(candi_bed6[,4])
candi_gtf=data.frame(0,0,0,0,0,0,0,0,0)
colnames(candi_gtf)=c("V1","V2","V3","V4","V5","V6","V7","V8","V9")
for(h in 1:length(um)){
subcandi_bed6=candi_bed6[which(candi_bed6[,4]==um[h]),]
strand=subcandi_bed6[1,6]
subcandi_bed6=peakreference(data = subcandi_bed6, merge = TRUE, overlap = 1)
subcandi_gtf=data.frame(subcandi_bed6[,1],"improve","exon",subcandi_bed6[,2]+1,subcandi_bed6[,3],"1000",strand,".",paste("gene_id ",um[h],"; transcript_id ",um[h],";",sep=""))
colnames(subcandi_gtf)=c("V1","V2","V3","V4","V5","V6","V7","V8","V9")
subcandi_mrna=data.frame(subcandi_bed6[1,1],"improve","transcript",min(subcandi_bed6[,2])+1,max(subcandi_bed6[,3]),"1000",strand,".",paste("gene_id ",um[h],";",sep=""))
colnames(subcandi_mrna)=c("V1","V2","V3","V4","V5","V6","V7","V8","V9")
candi_gtf=rbind(candi_gtf,subcandi_mrna,subcandi_gtf)
}
candi_gtf=candi_gtf[-1,]
write.table(candi_gtf,"JAAG.gtf",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
