##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep 26 14:51:15 2018
##########################################################################
##########################################################################
peak.files = list.files(path = "../../R6548_atac/Peaks/macs2",
                       pattern = "*.xls", full.names = TRUE)
peak.files = c(peak.files, list.files(path = "../../R6729_atac/Peaks/macs2", 
                        pattern = "*.xls", full.names = TRUE))
# import tbx peaks
peak.files = c(peak.files, "../data/tbx_90min_peaks_merged_macs2_p_5_filtered_N2_gene_assignment_TSS_WBcel235_analysis_v2.bed")

outDir = paste0("../results/DB_ABa_ABp")
if(!dir.exists(outDir)) dir.create(outDir);
manual.modifySampleInfos = TRUE

###############################
# libraries and functions
###############################
library("ChIPseeker");
library("rtracklayer")
#library(UpSetR);library("ChIPpeakAnno")
library("Vennerable") #install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
library("ggplot2");
library("GenomicFeatures");
library("GenomicRanges")
#library("DiffBind");library("IRanges");
source('functions_chipSeq.R')

###############################
# make design matrix and import peak files
###############################
source('functions_chipSeq.R')
design.matrix = make.design.matrix.from.peaks.files(peak.files = peak.files)

if(manual.modifySampleInfos){
  design.matrix$condition[which(design.matrix$condition=="60min")] = "90min"
  design.matrix$condition[which(design.matrix$condition=="140min")] = "200min"
  design.matrix$factor[grep("Aba", design.matrix$factor)] = "ABa"
  design.matrix$factor[grep("Abp", design.matrix$factor)] = "ABp"
  design.matrix$factor.condition = paste0(design.matrix$condition, "_", design.matrix$factor)
  design.matrix = design.matrix[order(design.matrix$condition, design.matrix$factor), ]
}

cat("-- import peak files as GRanges objects \n")
peaks.list = c()
for(n in 1:nrow(design.matrix)){
  #n = 14
  cat(n, '--', design.matrix$factor.condition[n], ":",  design.matrix$file.path[n],  '\n')
  xx = readPeakFile(design.matrix$file.path[n], as = "GRanges");
  
  if(seqlevelsStyle(xx) != "UCSC") seqlevelsStyle(xx) = "UCSC";
  
  peaks.list = c(peaks.list, xx)
  #eval(parse(text = paste0("pp.", n, " = xx")))
}

###############################
# peak overlapping checking
###############################
cat("-- compare peak overlapping and make plots \n")
#kk = c(1:6)

pdf(paste0(outDir, "/Comparison_peaks_for_ABa_ABp_overlapping_atac_tbx.pdf"),
    width = 12, height = 8)
source("functions_chipSeq.R") #
#par(mfrow=c(1,2))
Comparison.overlapping.peaks(design.matrix, peaks.list, toCompare="factor.condition", PLOT.p10 = TRUE)

dev.off()

########################################################
########################################################
# Section : define lineage specific peaks using called peaks
# define ABa- and ABp- specific peaks by pooling replicates and two time points (90 and 200 min)
########################################################
########################################################
library(rtracklayer)

sels = which(design.matrix$condition != "330min" & design.matrix$factor != 'tbx')
design.sels = design.matrix[sels,  ]
peaks.sels = peaks.list[sels]

conds = unique(design.sels$factor)
peaks = c()
peaknames = conds
pval = 10

for(fac in conds){
  # fac = 'ABa'
  kk = which(design.sels$factor == fac)
  peaks.merged = c() 
  
  for(k in kk) 
  {
    # k = 2
    xx = readPeakFile(design.sels$file.path[k], as = "GRanges");
    if(seqlevelsStyle(xx) != "UCSC") seqlevelsStyle(xx) = "UCSC";
    with.p.values = "X.log10.pvalue." %in% colnames(mcols(xx))
    p10 <- xx[mcols(xx)[,"X.log10.pvalue."] > pval]
    
    if(length(peaks.merged) == 0) {
      peaks.merged = reduce(p10, ignore.strand = TRUE)
    }else{
      peaks.merged = reduce(union(peaks.merged, p10, ignore.strand = TRUE), ignore.strand = TRUE)
    }
    cat(fac, ' -- all peaks :', length(xx), ' -- peaks <10-10: ', length(p10), ' -- merged peaks', length(peaks.merged), '\n')
  }
  peaks = c(peaks, peaks.merged)
}

pdf(paste0(outDir, "/Comparison_ATAC_peaks_for_early_ABa_ABp.pdf"),
    width = 12, height = 8)
source("functions_chipSeq.R") 

ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=c('early', 'ABa', 'ABp'), connectedPeaks="keepAll", main='early.vs.ABa.vs.ABp')

v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))

dev.off()

Define.Groups.peaks = FALSE
if(Define.Groups.peaks){
  
  ##########################################
  # define 7 groups of peaks using ABa ABp and 2to8cell.stage 
  ##########################################
  names(peaks) = conds
  p.all = union(peaks[[1]], peaks[[2]], ignore.strand = TRUE)
  p.all = union(p.all, peaks[[3]], ignore.strand = TRUE)
  p.all = reduce(p.all, ignore.strand = TRUE)
  
  export(p.all, con = paste0("../results/peakGroups/early_ABa_ABp_mergedPeaks.bed"))
  export(peaks[[1]], con = paste0("../results/peakGroups/early_allPeaks.bed"))
  export(peaks[[2]], con = paste0("../results/peakGroups/ABa_allPeaks.bed"))
  export(peaks[[3]], con = paste0("../results/peakGroups/ABp_allPeaks.bed"))
  
  
  p = p.all[overlapsAny(p.all, peaks[[1]])]
  p = p[overlapsAny(p, peaks[[2]])]
  p = p[overlapsAny(p, peaks[[3]])]
  
  export(p, con = paste0("../results/peakGroups/early_ABa_ABp_shared.bed"))
  
  p = p.all[!overlapsAny(p.all, peaks[[2]])]
  p = p[!overlapsAny(p, peaks[[3]])]
  
  export(p, con = paste0("../results/peakGroups/early_unique.bed"))
  
  p = p.all[overlapsAny(p.all, peaks[[2]])]
  p = p[overlapsAny(p, peaks[[1]])]
  p = p[!overlapsAny(p, peaks[[3]])]
  
  export(p, con = paste0("../results/peakGroups/early_ABa_shared.bed"))
  
  p = p.all[overlapsAny(p.all, peaks[[3]])]
  p = p[overlapsAny(p, peaks[[1]])]
  p = p[!overlapsAny(p, peaks[[2]])]
  
  export(p, con = paste0("../results/peakGroups/early_ABp.bed"))
  
  p = p.all[overlapsAny(p.all, peaks[[2]])]
  p = p[!overlapsAny(p, peaks[[1]])]
  p = p[!overlapsAny(p, peaks[[3]])]
  export(p, con = paste0("../results/peakGroups/ABa_unique.bed"))
  
  p = p.all[overlapsAny(p.all, peaks[[3]])]
  p = p[!overlapsAny(p, peaks[[2]])]
  p = p[!overlapsAny(p, peaks[[1]])]
  export(p, con = paste0("../results/peakGroups/ABp_unique.bed"))
  
  p = p.all[overlapsAny(p.all, peaks[[3]])]
  p = p[overlapsAny(p, peaks[[2]])]
  p = p[!overlapsAny(p, peaks[[1]])]
  export(p, con = paste0("../results/peakGroups/ABa_ABp_shared.bed"))
  
  ##########################################
  # define ABa and ABp specific peaks using ABa and ABp peaks
  ##########################################
  ABa = peaks.list[[which(design.matrix$factor.condition == "Aba_90min")]]
  ABp = peaks.list[[which(design.matrix$factor.condition == "Abp_90min")]]
  p = ABa 
  p <- p[mcols(p)[,"X.log10.pvalue."] > pval]
  p <- p[!overlapsAny(p, ABp)]
  
  export(p, con = paste0("results/motif_analysis/peaks_bed/ABa_unique_peaks.bed"))
  
  p = ABp
  p <- p[mcols(p)[,"X.log10.pvalue."] > pval]
  p <- p[!overlapsAny(p, ABa)]
  export(p, con = paste0("results/motif_analysis/peaks_bed/ABp_unique_peaks.bed"))
  ###############################
  # run motif analysis for those two unique peak sets 
  ###############################
  
}




