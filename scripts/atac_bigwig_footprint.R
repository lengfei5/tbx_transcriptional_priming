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

OutDir = "../data/bigWigs_PE/"
if(!dir.exists(OutDir)) dir.create(OutDir)

bamlist = list.files(path = "../data/Bams", pattern = "*.bam$", full.names = TRUE)

for(n in c(1:length(bamlist)))
{
  # n = 14
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
    if(!grepl('tbx', bam)){
      ga = readGAlignmentPairs(bam)
      #ga = readGAlignmentPairs(bam,param = ScanBamParam(flag=scanBamFlag(isDuplicate =FALSE))
    }else{
      ga = readGAlignments(bam)
    }
    
    ss = length(ga)/2
    xx = coverage(granges(ga))/s*10^6
    
    #xx = log2(xx+2^-6)
    export.bw(xx, con = paste0(OutDir, bw.name))
  }
}

########################################################
########################################################
# Section : split reads into nucleosome free and mono- dinucleosome regions
# initial code found 
# https://rockefelleruniversity.github.io/RU_ATAC_Workshop.html#greenleaf-dataset---finding-open-regions.
########################################################
########################################################
bamDir = "../data/Bams"
test.Nucleosome.Free = TRUE
Save.splitted.Bam = FALSE
Save.splitted.BW = FALSE
plot.fragmentSize = FALSE

outDir = "../results/bigWigs_cutSites_repAll_NFR"
if(!dir.exists(outDir)) dir.create(outDir)
outBam = paste0(outDir, "/bam_NFR")
outBW = paste0(outDir, "/bw_NFR")
if(!dir.exists(outBam)) dir.create(outBam)
if(!dir.exists(outBW)) dir.create(outBW)

bamlist = list.files(path = bamDir, pattern = "*.bam$", full.names = TRUE)
bamlist = bamlist[grep("tbx", bamlist, invert = TRUE)]

cat("directory to save -- ", outDir, "\n")

library(Rsubread)
library(GenomicAlignments)
library(rtracklayer)
library(magrittr)
library(dplyr)
library(ggplot2)

for(n in 1:length(bamlist))
{
  # n = 3
  cat("bam file -- ", bamlist[n], "\n")
  bam = bamlist[n]
  bname = paste0(gsub(".bam", "", basename(bam)))
  cat("bam name -- ", bname, "\n")
  
  # test chrV:10,643,153-10,649,406
  atacReads <- readGAlignmentPairs(bam,
                                   param = ScanBamParam(mapqFilter = 30,
                                                        flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isDuplicate =FALSE), 
                                                        what = c("qname", "mapq", "isize"), 
                                                        which = GRanges("chrV", IRanges(10000000, 12649406))))
  length(atacReads)
  atacReads
  ss = length(atacReads)/2
  
  if(test.Nucleosome.Free){
    nn = c(1:3)
    # Retrieving insert sizes
    atacReads_read1 <- GenomicAlignments::first(atacReads)
    insertSizes <- abs(elementMetadata(atacReads_read1)$isize)
    head(insertSizes)
    
    if(plot.fragmentSize){
      fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes, 
                                                                  Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)), 
                                                                                           Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) + 
        geom_line()
      
      fragLenPlot + theme_bw()
      fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
      fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") +
        geom_vline(xintercept = c(100), colour = "darkgreen") + 
        theme_bw()
      
    }
    # split the reads for NFR and mono
    atacReads_Open <- atacReads[insertSizes < 100, ]
    atacReads_MonoNuc <- atacReads[insertSizes > 180 & insertSizes < 240, ]
    atacReads_diNuc <- atacReads[insertSizes > 315 & insertSizes < 437, ]
    
    if(Save.splitted.Bam){
      # save splitted bam
      allRegionBam = paste0(outBam, "/", bname,  "_allRegions.bam")
      openRegionBam <- paste0(outBam, "/", bname, "_openRegions.bam")
      monoNucBam <- paste0(outBam, "/", bname, "_monoNuc.bam") 
      diNucBam <- paste0(outBam, "/", bname, "_diNuc.bam")
      
      export(atacReads, allRegionBam, format = "bam")
      export(atacReads_Open, openRegionBam, format = "bam")
      export(atacReads_MonoNuc, monoNucBam, format = "bam")
      export(atacReads_diNuc, diNucBam, format = "bam")
    }
    if(Save.splitted.BW){
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
  }else{
    nn = c(1)
  }
  ########################################################
  ########################################################
  # Section : Cutting sites from ATAC-seq data
  # ATAC-seq should generate shorter fragments (our nucleosome free regions) around smaller protected areas such as transcription factor binding sites.
  # We can therefore look for the pile-up of cut-sites around motifs of interest within different tissues/celltypes/samples.
  # To produce cut-sites from our BAM file we first resize our reads to 1bp and make the shift of 4/-5 bp depending on strand to 
  # adjust for expected shift from insertion of Tn5 transposase.
  # Here we will identify CTCF motifs passing an arbitary cut-off and then use soGGi to plot cut-sites around them
  ########################################################
  ########################################################
  #library(MotifDb)
  #library(Biostrings)
  #library(BSgenome.Hsapiens.UCSC.hg19)
  for(kk in nn){
    if(kk == 1){Reads = atacReads; bwName = paste0(outDir, '/', bname, ".bw") }
    if(kk == 2){Reads = atacReads_Open; bwName = paste0(outDir, '/', bname, "_Open.bw")}
    if(kk == 3){Reads = atacReads_MonoNuc; bwName = paste0(outDir, '/', bname, "_MonoNuc.bw")}
    cat(bwName, '\n')
    
    read1 <- first(Reads)
    read2 <- second(Reads)
    Firsts <- resize(granges(read1), fix = "start", 1)
    First_Pos_toCut <- shift(granges(Firsts[strand(read1) == "+"]), 4)
    First_Neg_toCut <- shift(granges(Firsts[strand(read1) == "-"]), -5)
    
    Seconds <- resize(granges(read2), fix = "start", 1)
    Second_Pos_toCut <- shift(granges(Seconds[strand(read2) == "+"]), 4)
    Second_Neg_toCut <- shift(granges(Seconds[strand(read2) == "-"]), -5)
    
    test_toCut <- c(First_Pos_toCut, First_Neg_toCut, Second_Pos_toCut, Second_Neg_toCut)
    cutsCoverage <- coverage(test_toCut)/ss*10^6
    export(cutsCoverage, con = bwName)
  
  }
  
}