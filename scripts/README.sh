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

# run rgt-hint for ABa and visualize the bias corrected signals
rgt-hint footprinting --atac-seq --paired-end --organism=ce11 --output-location=./results/footprint_hint --output-prefix=ABa_90min_chrV_test2 /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/data/test_footprint/ABa_90min_71327_chrV.bam /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/peakGroups/early_ABa_ABp_pooledPeaks_chrV.bed  

ml load libpng/1.2.58 # require old library lippng
rgt-hint tracks --bc --bigWig --organism=ce11 /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/data/test_footprint/ABa_90min_71327_chrV.bam /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/peakGroups/early_ABa_ABp_pooledPeaks_chrV.bed --output-location=./results/footprint_hint --output-prefix=ABa_90min_chrV_bw

## test single-end option
rgt-hint footprinting --atac-seq --organism=ce11 --output-location=./results/footprint_hint --output-prefix=ABa_90min_chrV_SEtest /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/data/test_footprint/ABa_90min_71327_chrV.bam /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/peakGroups/early_ABa_ABp_pooledPeaks_chrV.bed

# test wellington footprint
wellington_footprints.py /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/peakGroups/early_ABa_ABp_pooledPeaks_chrV_v2.bed /groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/data/test_footprint/ABa_90min_71327_chrV.bam ./results/footprint_wellington
    
