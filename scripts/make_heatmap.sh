##########
# bigwig and bed files (directories) required as inputs
# make heatmap and profiles using deeptools 
##########
cwd=`pwd`;
nb_cores=4;
OUT="$PWD/Results/heatmaps"
DIR_bw="$PWD/DATA/BigWigs_merged_chipseq_thomas_modEncode"
DIR_peak="$PWD/Results/peaks_targets_primed_direct"

#peak2use="10_10_PEAK_Union.bed"
peak="`ls ${DIR_peak}/*.bed |grep primed.bed` `ls ${DIR_peak}/*.bed |grep activated` `ls ${DIR_peak}/*.bed |grep undetectable` `ls ${DIR_peak}/*.bed |grep repressed` ";
ffs="`ls ${DIR_bw}/*.bw |grep tbx` `ls ${DIR_bw}/*.bw |grep -v tbx`"
cd $DIR_bw;
Samples=`ls *.bw|cut -d"_" -f1`
cd $pwd;

## parameters for deeptools
binsize=100;
dist=1000;

echo $peak;
#echo $Samples;
echo $ffs
#exit;

out0="HistoneMarks_TBX_within_chippeak_tbx"
#out_gcc=${out0}_peakExcluded_GCc
mkdir -p $OUT;
mkdir -p ${PWD}/logs;

cd $OUT
## make one big heatmap with defined grouped genomic regions in bed 
matrix_name="${out0}_matrix4heatmap.gz"
matrix2save="${out0}_matrix4save.txt"
heatmap="${out0}_heatmap_profile_RPKM"

qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N heatmap " module load python/2.7.3; module load deeptools/2.5.0.1-python2.7.3;\
if [ ! -e $matrix2save ]; then\
 computeMatrix reference-point --referencePoint=center -S "$ffs" -R $peak -a $dist -b $dist --outFileName "$matrix_name" --outFileNameMatrix "$matrix2save"\
 --skipZeros --binSize=$binsize  --sortUsingSamples 1 --sortUsing=mean; \
fi; \
plotHeatmap -m "$matrix_name" --colorMap RdBu\
  --sortRegions descend --sortUsing max -o ${heatmap}.pdf \
 --heatmapHeight 16\
 --refPointLabel peak --regionsLabel primed activated undetectable repressed; \
plotProfile -m "$matrix_name" -o ${heatmap}.png --plotHeight 12 --plotWidth 15 --colors green blue orange red; " 

cd $cwd 





 



