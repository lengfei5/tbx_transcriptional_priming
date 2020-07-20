#####################################
# this script is to do motifs analysis 
# using MEME suite
# the meme 4.12.0 is loaded from the module (current version)
#  
####################################
ml load bedtools/2.25.0-foss-2018b
ml load meme/5.0.4-foss-2018b-python-3.6.6
ml load xml-libxml/2.0132-gcccore-7.3.0-perl-5.28.0
#source activate meme
#conda env config vars set OMPI_MCA_opal_cuda_support=true
#source activate base
#module load meme/5.0.4-foss-2018b-python-2.7.15
#ml load python/3.6.6-foss-2018b
#ml load meme/5.0.4-foss-2018b-python-3.6.6
#MEME_path="/groups/cochella/jiwang/local/meme/bin/"

cwd=`pwd`;

# bed file as input
resDir='/groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/paper_revision/meme_res_test'
input_peaks="/groups/cochella/jiwang/Projects/Ariane/tbx_transcriptional_priming/results/paper_revision/peak_groups"

# output folder
out_fasta=${resDir}/seqs_fasta
out_ame=${resDir}/ame
out_meme=${resDir}/meme_chip

#nb_cores=2;
pwms="/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015.meme"
genome="/groups/bell/jiwang/Genomes/C_elegans/ce11/ce11_sequence/genome.fa"
bg="/groups/cochella/jiwang/Databases/motifs_TFs/background_cel_promoters/Promoter_sequences_all_ce11.fa"

primMotif="/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/tbx_motif.meme"
pval=0.0001;

mkdir -p $out_ame
mkdir -p $out_meme
mkdir -p $out_fasta

#mkdir -p $PWD/logs
#qsub -q public.q -o $cwd/logs -j yes -pe smp $nb_cores -cwd -b y -shell y -N fimo_scan "module load meme/4.9.1; fimo --bgfile $bg --oc $OUT --thresh $pval $pwms $sequence"
#module load meme/4.9.1;

for bb in `ls $input_peaks/*.bed `
do 
    #echo $bb;
    filename=`basename $bb`
    filename="${filename%.*}"
    #echo $filename
    if [ ! -e "${DIR}/${filename}.fa"  ]; then
	bedtools getfasta -fi $genome -bed $bb -name -fo ${out_fasta}/${filename}.fa 
    fi
done

for ff in `ls ${out_fasta}/*.fa`
do 
    filename=`basename $ff`
    filename="${filename%.*}"
    echo $filename
    
    meme-chip -db $pwms $ff -oc ${out_meme}/${filename} -fimo-skip -spamo-skip;   
    #ame --oc ${out_ame}/${filename} --control $bg --scoring totalhits --method fisher --bgformat 1 $ff $pwms;

    break
done
