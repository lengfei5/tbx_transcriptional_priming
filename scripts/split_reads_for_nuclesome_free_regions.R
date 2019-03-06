##########################################################################
##########################################################################
# Project:
# Script purpose: split the reads for nuclesome-free regsion (NFR) using pacakge "ATACseqQC" but it does not work well
# the current version is using the example from 
# https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Sep 24 15:45:57 2018
##########################################################################
##########################################################################
bamDir = "../../R6729_atac/alignments/BAMs_unique_rmdup"

outDir = "../results/splited_NFR"
if(!dir.exists(outDir)) dir.create(outDir)
outBam = paste0(outDir, "/bam_NFR")
outBW = paste0(outDir, "/bw_NFR")
if(!dir.exists(outBam)) dir.create(outBam)
if(!dir.exists(outBW)) dir.create(outBW)

bamlist = list.files(path = bamDir, pattern = "*.bam$", full.names = TRUE)
#bamlist = bamlist[grep("90min", bamlist, invert = TRUE)]

Save.splitted.Bam = TRUE;
save.splitted.BW = TRUE

library(Rsubread)
#pmapped <- propmapped(sortedBAM)
#pmapped
library(GenomicAlignments)
library(rtracklayer)
library(magrittr)
library(dplyr)
library(ggplot2)

for(n in 1:length(bamlist))
{
  # n = 1
  cat("bam file -- ", bamlist[n], "\n")
  bam = bamlist[n]
  
  bname = paste0(gsub(".bam", "", basename(bam)))
  cat("bam name -- ", bname, "\n")
  
  #splitDir = paste0("results/splited_NFR/", bname)
  
  cat("directory to save -- ", outDir, "\n")
  
  sortedBAM = bam
  
  atacReads <- readGAlignmentPairs(sortedBAM, 
                                   param = ScanBamParam(mapqFilter = 30, 
                                                        flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isDuplicate =FALSE), 
                                                        what = c("qname", "mapq", "isize")))
  # length(atacReads)
  #atacReads
  
  # Retrieving insert sizes
  atacReads_read1 <- GenomicAlignments::first(atacReads)
  insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
  head(insertSizes)
  
  fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
                                                              Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
                                                                                       Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
    geom_line()
  
  fragLenPlot + theme_bw()
  fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
  fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") +
    geom_vline(xintercept = c(100), colour = "darkgreen") + 
    theme_bw()
  
  # split the reads for NFR and mono
  atacReads_Open <- atacReads[insertSizes < 100, ]
  atacReads_MonoNuc <- atacReads[insertSizes > 180 & insertSizes < 240, ]
  atacReads_diNuc <- atacReads[insertSizes > 315 & insertSizes < 437, ]
  
  # save splitted bam
  allRegionBam = paste0(outBam, "/", bname,  "_allRegions.bam")
  openRegionBam <- paste0(outBam, "/", bname, "_openRegions.bam")
  monoNucBam <- paste0(outBam, "/", bname, "_monoNuc.bam") 
  diNucBam <- paste0(outBam, "/", bname, "_diNuc.bam")
  
  export(atacReads, allRegionBam, format = "bam")
  export(atacReads_Open, openRegionBam, format = "bam")
  export(atacReads_MonoNuc, monoNucBam, format = "bam")
  export(atacReads_diNuc, diNucBam, format = "bam")
  
  # save splitted bigwig
  allRegionBW = paste0(outBW, "/", bname, "_allRegions.bw")
  openRegionBW <- paste0(outBW, "/", bname, "_openRegions.bw")
  monoNucBW <- paste0(outBW, "/", bname, "_monoNuc.bw") 
  diNucBW <- paste0(outBW, "/", bname, "_diNuc.bw")
  
  #coverage(granges(ga))
  export.bw(coverage(granges(atacReads)), allRegionBW)
  export.bw(coverage(granges(atacReads_Open)), openRegionBW)
  export.bw(coverage(granges(atacReads_MonoNuc)), monoNucBW)
  export.bw(coverage(granges(atacReads_diNuc)), diNucBW)
  
}

########################################################
########################################################
# Section: try ATACseqQC R pacakge which does not work at all
########################################################
########################################################
Use.ATACseqQC = FALSE
if(Use.ATACseqQC){
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
}