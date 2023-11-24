# Join and judge Analysis for Annotating coding Genes (JAAG)
## Summarize
## Requirements
minimap2 >= 2.16

samtools >= 1.6

python >= 2.7.10

diamond >= 2.0.4

BEDTools >= 2.27.1

R >= 3.6.0
## Getting start
### 0. Data Preparation (input)
ref.fa : genome reference sequences (.fa)

movieX.flnc.fa: Full Length Non-Chimeric (FLNC) read file (.fa)

### 1. ISO-seq-based
#### i. Aligning Iso-Seq data
```
minimap2 -x map-pb -d ref.mmi ref.fa
minimap2 -t 8 -ax splice -uf --secondary=no -C5 -a  refgen.mmi movieX.flnc.fa >  movieX.sam
samtools sort movieX.sam > movieX.sort.bam
samtools view -h -o movieX.sort.sam movieX.sort.bam
```
#### ii. FLNC reads collapse, merge, and filtering
##### Collapse multiple transcriptional reads into single transcript model
```
python tama/tama_collapse.py -s movieX.sort.sam -f ref.fa -p prefix1 -x capped -a 100 -z 100 -sj sj_priority -sjt 20 -lde 2 -log log_off
```
-x [capped/ no_cap]: Iso-Seq library preparation is performed with 5' cap selection or not;

-a [integer] (default 10): The amount of tolerance at the 5' end of the reads for grouping together.

-z [integer] (default 10): The amount of tolerance at the 3' end of the reads for grouping together.

-sj [no_priority/ sj_priority]: The parameter decides if use splice junction priority information for final model.

-sjt [integer] (default 10): The length threshold for the regions to observe for splice junction priority.

-lde [integer] (default 1000): The number of errors that is allowed on each side of the splice junction within the length specified by the -sjt argument.

##### Merge multiple transcriptomes into single annotation (If Iso-Seq data has multiple runs)
```
python tama/tama_merge.py -f merge.txt -p merged_annos
```
-p Output prefix

-f flielist txt file (created by users)

The example of filelist file
```
prefix_1.bed  capped  1,1,1  prefix_1
prefix_2.bed  capped  1,1,1  prefix_2
```
This format includes 4 column without header: column 1 is the output bed file from tama collapse step; column 2 indicates if this library is capped or not; column 3 ranks the confidence of the information from each library, from the first to the last, they are start site, splice junctions, and end sites; column 4 is used for prefix name when output.

##### Counting read support and filtering
###### Count read support information from each TAMA Collaspe run
```
python tama/tama_read_support_levels.py -f read_support_filelist.txt -o merge -m merge.txt
```
-o Output prefix

-f flielist txt file (created by users)

-m merge file (Here, it refers to merge.txt file from TAMA Merge step)

The example of filelist file
```
prefix_1  prefix_1_trans_read.bed  trans_read
prefix_2  prefix_2_trans_read.bed  trans_read
```
This format includes 3 column without header: column 1 is the source name defined by user that can trace reads back to which library when merge together ; column 2 is the "trans_read.bed" file generated in TAMA Collaspe step; column 3 indicates the type of countting (In this step it is "trans_read"). 

###### Remove gene model with genomic polyA stretches at the end of three primer
```
python tama/tama_remove_polya_models_levels.py -b merged_annos.bed -f polya_filelist.txt -r merge_read_support.txt -o prefix -k keep_multi
```
-b Annotation bed file generated from merge step

-f Filelist file with Poly-A file names (created by users)

-r Read support file obtained from the read support step above

-o Output prefix

-k [keep_multi/ remove_multi] All multi-exon gene will not be removed in this step if keep_multi

The example of the polyA Filelist file
```
prefix_1  prefix_1_polya.txt
prefix_2  prefix_2_polya.txt
```
This format includes 2 column without header: column 1 is the source name defined by user that can trace reads back to which library when merge together ; column 2 is the "polya.txt" file generated in TAMA Collaspe step.

###### Recount read support information after polyA filtering
```
python tama/tama_read_support_levels.py -f read_support_rm_polya.txt -o rm_polya -m merge_rmpolya_polya_report.txt -mt filter
```
-f flielist txt file (created by users)

-o Output file prefix

-m merge file (Here, it refers to merge_rmpolya_polya_report.txt file from polyA-transcripts removal step)

-mt [tama cupcake filter] "filter" indicates the input file is produced from one of TAMA filtering steps.

The example of filelist file (?????)
```
prefix_1  prefix_1_trans_read.bed  trans_read
prefix_2  prefix_2_trans_read.bed  trans_read
```

###### Remove transcript with only single read support
```
python tama/tama_remove_single_read_models_levels.py -b merge_rmpolya.bed -r rm_polya_read_support.txt -l transcript -k keep_multi -o merge_rmpolya_rmsingle
```
-b Annotation bed file which is obtained from polyA filtering step
-r Read support file from recount read support step after polyA filtering
-o Output prefix (required)
-l [gene transcript] Level of removal
-k [keep_multi remove_multi] Keep or remove multi-exon models if the transcript with only one read support.


#from tama-based.bed12 to tama-based.coding.bed12
module load igmm/apps/miniconda3/4.5.11
source activate /home/s1874451/schoenebeck_group/WENGANG/anaconda_env/diamond
diamond makedb --in input_db.fa -d input_db
diamond blastx --db input_db --query merge_rmpolya_rmsingle.fa  --sensitive -f 6 --out merge_rmpolya_rmsingle.input_db.results --max-hsps 0 -c1 --strand plus
module load igmm/apps/BEDTools/2.27.1
bed12ToBed6 -i merged_rmpolya_rmsingle.bed > merged_rmpolya_rmsingle.bed6
run JAAG_readthrough.V3.r (protein-guided) to get transcripts rename clue file (merged_rmpolya_rmsingle_coding.readthrough.reassign) and rename bed file (merged_rmpolya_rmsingle_coding.readthrough_reassign.bed)
```


### 2. RNA-seq-based
### 3. JAAG run
