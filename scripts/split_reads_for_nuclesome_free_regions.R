##########################################################################
##########################################################################
# Project:
# Script purpose: split the reads for nuclesome-free regsion (NFR) using pacakge "ATACseqQC"
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Sep 24 15:45:57 2018
##########################################################################
##########################################################################
install.pkgs = FALSE
if(install.pkgs)
{
  source("https://bioconductor.org/biocLite.R")
  biocLite("BSgenome.Celegans.UCSC.ce11")
  biocLite("TxDb.Celegans.UCSC.ce11.refGene")
  biocLite("Rsubread")
  biocLite("Rsamtools")
  biocLite("ATACseqQC")
}

library("ATACseqQC")
library("Rsubread")
library("Rsamtools")
library(GenomicAlignments)
#library("ggplot2")
#library("devtools")
#library("magrittr")
library(BSgenome.Celegans.UCSC.ce11)
library(TxDb.Celegans.UCSC.ce11.refGene)

bamlist = list.files(path = "../R6548_atac/alignments/BAMs_unique_rmdup", pattern = "*.bam$", full.names = TRUE)

txs <- transcripts(TxDb.Celegans.UCSC.ce11.refGene)
#txs <- txs[seqnames(txs) %in% "chrI"]

genome <- Celegans

#dir.create(outPath)

for(n in c(1))
{
  # n = 1
  cat("bam file -- ", bamlist[n], "\n")
  bam = bamlist[n]
  
  #ga = readGAlignmentPairs(bam)
  #ga = readGAlignmentPairs(bam,param = ScanBamParam(flag=scanBamFlag(isDuplicate =FALSE))
  #xx = coverage(granges(ga))
  #export.bw(xx, con = paste0(, bwname))

  bname = paste0(gsub(".bam", "", basename(bam)))
  cat("bam name -- ", bname, "\n")
  
  splitDir = paste0("results/splited_NFR/", bname)
  if(!dir.exists(splitDir)) dir.create(splitDir)
  
  cat("directory to save -- ", splitDir, "\n")
  
  #seqlev <- "chrI" ## subsample data for quick run
  #which <- as(seqinfo(Celegans)[seqlev], "GRanges")
  #tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")
  
  bb = readBamFile(bam, asMates = TRUE)
  
  cat("start to shift the reads ----\n")
  
  bb <- shiftGAlignmentsList(bb)
  
  #bb.shifted <- file.path(splitDir,  "_shifted.bam"))
  
  # export(bb, bb.shifted)
  #names(bb) <- mcols(bb)$qname
  
  cat("start to split reads into nucleosome free and ... \n")
  
  objs <- splitGAlignmentsByCut(bb, txs=txs, genome=genome)
  
  ### Save the binned alignments into bam files
  cat("start to save the bam files ...\n")
  null <- writeListOfGAlignments(objs, splitDir)
  
}






