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
### 0. Data Preparation
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
```
python tama/tama_collapse.py -s movieX.sort.sam -f ref.fa -p prefix -x capped -a 100 -z 100 -d merge_dup -sj sj_priority -sjt 20 -lde 2 -log log_off
python tama/tama_merge.py -f merge.txt -p merged_annos
#TAMA Read Support (https://github.com/GenomeRIK/tama/wiki/TAMA-GO:-Read-Support)   see TXT.format
python ~/WGsoftware/tama/tama_go/read_support/tama_read_support_levels.py -f read_support_filelist.txt -o merge -m _merge.txt
#TAMA filter_out (https://github.com/GenomeRIK/tama/wiki/TAMA-GO:-Transcript-Filtering)   see TXT.format
python ~/WGsoftware/tama/tama_go/filter_transcript_models/tama_remove_polya_models_levels.py -b merged_annos.bed -f polya_filelist.txt -r merge_read_support.txt -o prefix -k keep_multi
#TAMA Read Support for polyA-removed transcripts
python ~/WGsoftware/tama/tama_go/read_support/tama_read_support_levels.py -f read_support_rm_polya.txt -o rm_polya -m merge_rmpolya_polya_report.txt -mt filter
#TAMA Transcript Filtering (reads_support)
python ~/WGsoftware/tama/tama_go/filter_transcript_models/tama_remove_single_read_models_levels.py -b merge_rmpolya.bed -r rm_polya_read_support.txt -l transcript -k keep_multi -o merge_rmpolya_rmsingle

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
