# test tbx motif scannign around lsy-6
ml load bedtools/2.27.1-foss-2017a
cd /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/peakGroups
bedtools getfasta -fi /groups/cochella/jiwang/Genomes/C_elegans/ce11/ce11_sequence/genome.fa -bed testRegions_aroundlsy6.bed > ../motif_analysis/testRegions_aroundlsy6.fa

## add +1 to start coord. of the header to make it 1-based
cd /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/motif_analysis
awk -F ":" '{OFS=""; split($2,a,"-"); if(a[1]) print $1,":",a[1]+1,"-",a[2]; else print; }' testRegions_aroundlsy6.fa > testRegions_aroundlsy6_v1.fa

cd /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming
ml load meme/4.12.0-foss-2017a
# scan only tbx
fimo --oc results/motif_analysis/testRegions_tbx --parse-genomic-coord -thresh 0.001 /groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/tbx_motif.meme results/motif_analysis/testRegions_aroundlsy6_v1.fa
# scan all motif
fimo --oc results/motif_analysis/testRegions_all --parse-genomic-coord -thresh 0.0001 /groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015.meme results/motif_analysis/testRegions_aroundlsy6_v1.fa

# test rgt-hint footprint 
cd /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/data/test_footprint/HINT_ATACTest
conda activate RGT
rgt-hint footprinting --atac-seq ATAC.bam ATACPeaks.bed

# run rgt-hint for ABa
rgt-hint footprinting --atac-seq --paired-end --organism=ce11 --output-location=./results/footprint_hint --output-prefix=ABa_90min /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/data/Bams/ABa_90min_71327.bam /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/peakGroups/early_ABa_ABp_pooledPeaks_chrV.bed  
