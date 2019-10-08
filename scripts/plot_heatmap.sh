#!/bin/bash
#$ -q public.q
#$ -cwd
#$ -m ea
#$ -o $PWD/logs
#$ -e $PWD/logs
#$ -pe smp 4-10

#####################
# the script is to make figure 1A the panel for Rybp in WT
# for Jorge's manuscript
#####################
out_dir="/groups/bell/jiwang/Projects/Jorge/figures/FigS1E"

#bw_file="/groups/bell/jiwang/Projects/Jorge/NGS_requests/R6118_chipseq/bigWigs_GCc_merged/Cbx7_AN312_65721.65722.merged.bw"
bw_file="/groups/bell/jiwang/Projects/Jorge/NGS_requests/R6820_chipseq/bigWigs_GCc_merged/Rybp_WT.AN312_75656.75659.merged_mq_30.bw"
bw_file_2="Rybp_EedKO.515D10_75657.75660.merged_mq_30.bw"
#bw_file="/groups/bell/jiwang/Projects/Jorge/NGS_requests/R6820_chipseq/bigWigs_GCc_merged/Rybp_EedKO.515D10_75657.75660.merged_mq_30.bw"
#bw_file="/groups/bell/jiwang/Projects/Jorge/NGS_requests/R6820_chipseq/bigWigs_GCc_merged/Rybp_RybpKO.AN3E9_75658.75661.merged_mq_30.bw"
bed_file_1="/groups/bell/jiwang/Projects/Jorge/figures/All_Refseq_TSS_Analysis/sharedTSS_backup.bed"
bed_file_2="/groups/bell/jiwang/Projects/Jorge/figures/All_Refseq_TSS_Analysis/variantTSS_backup.bed"
#bed_file_3="/groups/bell/jiwang/Projects/Jorge/figures/All_Refseq_TSS_Analysis/noBioCapTSS.bed"
#in_dir2="/groups/bell/jiwang/Projects/Jorge/NGS_requests/R6118_chipseq/bigWigs_GCc_merged"
#file_name=`basename ${bw_file}`;
file_name="Rybp_WT_vs_EedKO"

mkdir -p $out_dir
module unload python/2.7.3;
#module unload pythonlib;
module unload deeptools;
module load python/2.7.3;
#module load pysam/0.10.0;
#module load deeptools/2.2.3-python2.7.3;
#module load pythonlib;
#module load deeptools/2.5.0.1-python2.7.3;
module load deeptools/2.2.3-python2.7.3
## only for Rybp in WT (AN312 cells)
computeMatrix reference-point \
-S $bw_file \
-R $bed_file_1 $bed_file_2 $bed_file_3 \
--referencePoint center \
--beforeRegionStartLength 5000 \
--afterRegionStartLength 5000 \
--sortRegions no \
-o ${out_dir}/WT_PcG-prots-Biocap_Refseq-Shared-Variant-TSS.gz


module load deeptools/2.5.0.1-python2.7.3
plotHeatmap -m ${out_dir}/WT_PcG-prots-Biocap_Refseq-Shared-Variant-TSS.gz \
--sortRegions no \
--zMin 0 --zMax 12 \
--regionsLabel "Shared" "Variant" "nonBioCap" \
--samplesLabel $file_name \
--plotTitle $file_name \
--colorMap YlGnBu \
--legendLocation lower-center \
-o ${out_dir}/${file_name}_Shared-Variant-TSS_heatmap.pdf

 plotProfile -m ${out_dir}/WT_PcG-prots-Biocap_Refseq-Shared-Variant-TSS.gz \
 	-out ${out_dir}/${file_name}_Shared-Variant-TSS_profile.pdf \
 	--colors red blue green \
  --regionsLabel "Shared" "Variant" "nonBioCap" \
	--samplesLabel $file_name \
 	--plotTitle $file_name \
  --yMin 0 --yMax 18\
  --numPlotsPerRow 3 \
  --legendLocation upper-right \
  --outFileNameData ${out_dir}/${file_name}.txt
