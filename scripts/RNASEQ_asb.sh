mkdir sambam
mkdir temp
mkdir gtf
for sample_name  in `cat rna_seq.list`; do
hisat2 -q -x ./g_index/ref \
 -1 ./rnaseq/${sample_name}_1.clean.fq.gz -2 ./rnaseq/${sample_name}_2.clean.fq.gz \
 -S ./sambam/${sample_name}.sam 2> ./sambam/${sample_name}.hisat2.log
samtools sort -o ./sambam/${sample_name}.sort.bam ./sambam/${sample_name}.sam
samtools index ./sambam/${sample_name}.sort.bam
rm ${sample_name}.sam
samtools view -H sambam/${sample_name}.sort.bam |grep "@SQ" |awk '{print $2}' |awk -F':' '{print $2}' > temp/${sample_name}.chr
for chr  in `cat ./temp/${sample_name}.chr`; do
samtools view -b sambam/${sample_name}.sort.bam ${chr} > sambam/${sample_name}.${chr}.sort.bam
samtools  view -L ./isoseq/merged_rmpolya_rmsingle_coding.readthrough_reassign.N.bed6 ./sambam/${sample_name}.${chr}.sort.bam |awk '$2!=83 && $2!=163 && $2!=161 && $2!=81 && $2!=419 && $2!=339 '|cut -f1 |sort |uniq > ./temp/${sample_name}.${chr}.sort.N.name
samtools  view -L ./isoseq/merged_rmpolya_rmsingle_coding.readthrough_reassign.P.bed6 ./sambam/${sample_name}.${chr}.sort.bam |awk '$2!=99 && $2!=147 && $2!=403 && $2!=355 && $2!=145 && $2!=97 '|cut -f1 |sort |uniq > ./temp/${sample_name}.${chr}.sort.P.name
samtools view -f256  ./sambam/${sample_name}.${chr}.sort.bam |cut -f1 |sort |uniq > ./temp/${sample_name}.${chr}.Multi.name
cat ./temp/${sample_name}.${chr}.sort.N.name ./temp/${sample_name}.${chr}.sort.P.name |sort|uniq > ./temp/${sample_name}.${chr}.name
awk -f ./scripts/vlookup.awk ./temp/${sample_name}.${chr}.Multi.name ./temp/${sample_name}.${chr}.name |grep "KP" | awk '{print $1}' > ./temp/${sample_name}.${chr}.RM.name
java -jar $PICARD/picard.jar FilterSamReads  \
I=./sambam/${sample_name}.${chr}.sort.bam  O=./sambam/${sample_name}.${chr}.JAAG.bam READ_LIST_FILE=./temp/${sample_name}.${chr}.RM.name  FILTER=excludeReadList
done
samtools merge sambam/${sample_name}.JAAG.bam sambam/${sample_name}*.JAAG.bam
rm sambam/${sample_name}.*.sort.bam sambam/${sample_name}.*.JAAG.bam sambam/${sample_name}.sort.bam.bai
stringtie -o ./gtf/${sample_name}.gtf  sambam/${sample_name}.JAAG.bam
#build bed file based on strand
awk '{if($7==".") print $0}' ./gtf/${sample_name}.gtf | awk '{print $12}' |sort |uniq  > ./gtf/${sample_name}.RM.SS.gtf.name
awk '{if($3=="exon"&&$7=="+") print $0}' ./gtf/${sample_name}.gtf |awk  'BEGIN{OFS="\t"} {print $1, $4, $5, $12}' > ./gtf/${sample_name}.P.bed6
awk '{if($3=="exon"&&$7=="-") print $0}' ./gtf/${sample_name}.gtf |awk  'BEGIN{OFS="\t"} {print $1, $4, $5, $12}' > ./gtf/${sample_name}.N.bed6
bedtools intersect -a ./gtf/${sample_name}.P.bed6 -b ./isoseq/merged_rmpolya_rmsingle_coding.readthrough_reassign.P.bed6 > ./gtf/${sample_name}.overlap.P.bed6
bedtools intersect -a ./gtf/${sample_name}.N.bed6 -b ./isoseq/merged_rmpolya_rmsingle_coding.readthrough_reassign.N.bed6 > ./gtf/${sample_name}.overlap.N.bed6
cut -f4 ./gtf/${sample_name}.overlap.N.bed6 |sort |uniq > ./gtf/${sample_name}.RM.N.gtf.name
cut -f4 ./gtf/${sample_name}.overlap.P.bed6 |sort |uniq > ./gtf/${sample_name}.RM.P.gtf.name
cat ./gtf/${sample_name}.RM.N.gtf.name ./gtf/${sample_name}.RM.P.gtf.name ./gtf/${sample_name}.RM.SS.gtf.name > ./gtf/${sample_name}.RM.gtf.name
awk -f ./scripts/vlookup2.awk  ./gtf/${sample_name}.RM.gtf.name ./gtf/${sample_name}.gtf |grep 'KP'|awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > ./gtf/${sample_name}.JAAG.gtf
rm ./gtf/*.bed6 ./gtf/*name temp/${sample_name}.* 
done
l gtf/*.JAAG.gtf |awk -F' ' '{print $9}' > gtf/gtf.list
stringtie --merge -o gtf/merge.candi.gtf gtf/gtf.list
