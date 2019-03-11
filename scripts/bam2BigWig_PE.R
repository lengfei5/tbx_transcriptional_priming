##########################################################################
##########################################################################
# Project:
# Script purpose: make bigwig file for paired_end bam
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep 26 14:22:26 2018
##########################################################################
##########################################################################
library(GenomicAlignments)
library(rtracklayer)

OutDir = "../data/bigWigs_PE_log2/"
if(!dir.exists(OutDir)) dir.create(OutDir)

bamlist = list.files(path = "../data/Bams", pattern = "*.bam$", full.names = TRUE)

for(n in c(1:length(bamlist)))
{
  # n = 1
  bam = bamlist[n]
  bw.name = basename(bam)
  bw.name = gsub(".bam", ".bw", bw.name)
  bw.name = gsub("_uniq_rmdup", '', bw.name)
  bw.name = gsub("140min", '200min', bw.name)
  bw.name = gsub("60min", '90min', bw.name)
  bw.name = gsub("Aba", 'ABa', bw.name)
  bw.name = gsub("Abp", 'ABp', bw.name)
  
  cat("bam file: ", bamlist[n], '-- ', "bw name: ", bw.name, "\n")

  if(! file.exists(paste0(OutDir, bw.name))){
    ga = readGAlignmentPairs(bam)
    #ga = readGAlignmentPairs(bam,param = ScanBamParam(flag=scanBamFlag(isDuplicate =FALSE))
    xx = coverage(granges(ga))
    xx = log2(xx+2^-6)
    export.bw(xx, con = paste0(OutDir, bw.name))
  }
}
