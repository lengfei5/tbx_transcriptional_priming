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

OutDir = "bigWigs_PE/"
if(!dir.exists(OutDir)) dir.create(OutDir)

bamlist = list.files(path = "../R6548_atac/alignments/BAMs_unique_rmdup", pattern = "*.bam$", full.names = TRUE)

for(n in c(1:length(bamlist)))
{
  # n = 1
  cat("bam file -- ", bamlist[n], "\n")
  bam = bamlist[n]
  bw.name = basename(bam)
  bw.name = gsub(".bam", ".bw", bw.name)
  ga = readGAlignmentPairs(bam)
  #ga = readGAlignmentPairs(bam,param = ScanBamParam(flag=scanBamFlag(isDuplicate =FALSE))
  xx = coverage(granges(ga))
  
  export.bw(xx, con = paste0(OutDir, bw.name))
  
}
