#############
# bigwig and bed files (directories) required as inputs
# make heatmap and profiles using deeptools 
##############
# all abosute paths here
OUT="/groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/heatmaps"
#DIR_bw="/groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/bigWigs_PE_toMerge"
DIR_bw='/groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/data/bigWigs_PE_log2'
DIR_peak="/groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/peakGroups"
out0="atacSeq_tbx_peakGroups"

cwd=`pwd`

ffs="${DIR_bw}/sorted.2to8cell.stage_0min_73805.bw \
${DIR_bw}/ABa_90min_71327.bw ${DIR_bw}/ABa_90min_73807.bw \
${DIR_bw}/ABp_90min_71328.bw ${DIR_bw}/ABp_90min_73806.bw \
${DIR_bw}/ABa_200min_71329.bw ${DIR_bw}/ABa_200min_73809.bw \
${DIR_bw}/ABp_200min_71330.bw ${DIR_bw}/ABp_200min_73808.bw \
${DIR_bw}/ABa_330min_71331.bw ${DIR_bw}/ABa_330min_73811.bw \
${DIR_bw}/ABp_330min_71332.bw ${DIR_bw}/ABp_330min_73810.bw"

peak="${DIR_peak}/early_ABa_ABp_pooledPeaks.bed"


#`ls ${DIR_bw}/*.bw|grep tbx_merged` "
#cd $DIR_bw;

## parameters for deeptools
binsize=10;
dist=500;
#echo 'peaks files : '
echo $peak;
#echo 'bw files : '
echo $ffs

#out_gcc=${out0}_peakExcluded_GCc
mkdir -p $OUT;
mkdir -p ${OUT}/logs;

## make one big heatmap with defined grouped genomic regions in bed 
matrix_name="${out0}_matrix4heatmap.gz"
matrix2save="${out0}_matrix4save.txt"
heatmap="${out0}_heatmap_profile_log2RPKM_all_kmean10"
#exit
cd $OUT

script=${out0}_deeptools.sh

# make submission script
cat <<EOF > $script
#!/usr/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --time=120
#SBATCH --mem=4000

#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH -o $OUT/logs/${out0}.out
#SBATCH -e $PWD/logs/${out0}.err
#SBATCH --job-name heatmap

module load python/2.7.13-foss-2017a
module load deeptools/3.0.1-foss-2017a-python-2.7.13

computeMatrix reference-point --referencePoint=center \
-S $ffs -R $peak \
-a $dist -b $dist \
--outFileName "$matrix_name" \
--outFileNameMatrix "$matrix2save" \
--skipZeros --binSize=$binsize --sortUsingSamples 1 --sortUsing=mean

plotHeatmap -m "$matrix_name" --colorMap RdBu \
--sortRegions descend --sortUsing max -o ${heatmap}.pdf \
--heatmapHeight 16 \
--heatmapWidth 2 \
--refPointLabel peak \
--kmeans 10

# plotProfile -m "$matrix_name" -o ${heatmap}.png --plotHeight 12 --plotWidth 15 --colors green blue orange red

EOF

cat $script;
sbatch $script

cd $cwd
