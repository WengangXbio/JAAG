library(stringr)
library(TCseq)
options("scipen"=15)
args <- commandArgs(trailingOnly = TRUE)
bed6file=read.table(args[1],head=F)
n=length(unique(bed6file[,4]))
genegtf=data.frame(0,0,0,0)
colnames(genegtf)=c("chr","start","end","id")
strand=c()
for(i in 1:n){
bed6file_gene=bed6file[which(bed6file[,4]==unique(bed6file[,4])[i]),]
merged_gene <- peakreference(data = bed6file_gene, merge = TRUE, overlap = 1)
merged_gene[,4]=bed6file_gene[1,4]
strand=c(strand,as.character(bed6file_gene[1,6]))
genegtf=rbind(genegtf,merged_gene)
}
genegtf=genegtf[-1,]
fwe_bed=data.frame(0,0,0,0,0,0,0,0,0,0,0,0)
colnames(fwe_bed)=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12")
for(j in 1:n){
genegtfsub=genegtf[which(genegtf[,4]==unique(genegtf[,4])[j]),]
a=genegtfsub[1,3]-genegtfsub[1,2]
b=0
if(dim(genegtfsub)[1]!=1)
for(e in 2:dim(genegtfsub)[1]){
a=paste(a,",",genegtfsub[e,3]-genegtfsub[e,2],sep="")
b=paste(b,",",genegtfsub[e,2]-genegtfsub[1,2],sep="")
}
uey=data.frame(V1=genegtfsub[1,1],V2=min(genegtfsub[,2]),V3=max(genegtfsub[,3]),V4=paste(genegtfsub[1,4],"-T",sep=""),V5=40,V6=strand[j],V7=min(genegtfsub[,2]),V8=max(genegtfsub[,3]),V9="255,0,0",V10=dim(genegtfsub)[1],V11=a,V12=b)
fwe_bed=rbind(fwe_bed,uey)
}
fwe_bed=fwe_bed[-1,]
write.table(fwe_bed,args[12],quote=FALSE,sep="\t",row.names=F,col.names=F)
