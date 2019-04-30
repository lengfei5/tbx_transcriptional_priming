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
design.matrix = make.design.matrix.from.file.list(peak.files = peak.files)

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
peaks.all = NULL
for(n in 1:nrow(design.matrix)){
  #n = 14
  cat(n, '--', design.matrix$factor.condition[n], ":",  design.matrix$file.path[n],  '\n')
  xx = readPeakFile(design.matrix$file.path[n], as = "GRanges");
  
  if(seqlevelsStyle(xx) != "UCSC") seqlevelsStyle(xx) = "UCSC";
  
  peaks.list = c(peaks.list, xx)
  if(n==1) {
    peaks.all = xx; 
  }else{
    peaks.all = union(peaks.all, xx, ignore.strand = TRUE)
  }
}

peaks.all = reduce(peaks.all)

export(peaks.all, con = paste0("../results/peakGroups/early_ABa_ABp_pooledPeaks.bed"))

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

#sels = which(design.matrix$condition != "330min" & design.matrix$factor != 'tbx' & design.matrix$condition != "200min")
sels = which(design.matrix$factor.condition == '90min_ABa' | design.matrix$factor.condition == '90min_ABp')
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
      peaks.merged = reduce(intersect(peaks.merged, p10, ignore.strand = TRUE), ignore.strand = TRUE)
    }
    cat(fac, ' -- all peaks :', length(xx), ' -- peaks <10-10: ', length(p10), ' -- merged peaks', length(peaks.merged), '\n')
  }
  peaks = c(peaks, peaks.merged)
}

pdf(paste0(outDir, "/Comparison_ATAC_peaks_for_ABa_ABp_replicateOverlapPeaks.pdf"),
    width = 12, height = 8)
source("functions_chipSeq.R") 

ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=c('ABa', 'ABp'), connectedPeaks="keepAll", main='ABa.vs.ABp')

v <- venn_cnt2venn(ol.peaks$vennCounts)
try(plot(v))

dev.off()

##########################################
# define 7 groups of peaks using ABa ABp and 2to8cell.stage 
##########################################
Define.Groups.peaks = FALSE
if(Define.Groups.peaks){
  names(peaks) = conds
  p.all = union(peaks[[1]], peaks[[2]], ignore.strand = TRUE)
  p.all = union(p.all, peaks[[3]], ignore.strand = TRUE)
  p.all = reduce(p.all, ignore.strand = TRUE)
  
  export(p.all, con = paste0("../results/peakGroups/early_ABa_ABp_pooledPeaks.bed"))
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

########################################################
########################################################
# Section: quantify signals within peaks and background regions
# Differential Binding analysis
########################################################
########################################################
DIR.bams = "../data/Bams"
DIR.peaks = "../results/peakGroups"

resDir = '../results/DB_ABa_ABp/'
tableDir = paste0(resDir, "tables/")
version.analysis = "_20190312"

addBackground = TRUE
run.featureCount.quantification = FALSE
run.peak2geneAssignment = FALSE

batch.removal = FALSE

if(!dir.exists(resDir)) system(paste0('mkdir -p ', resDir))
if(!dir.exists(tableDir)) system(paste0('mkdir -p ', tableDir))

bam.files = list.files(path = DIR.bams, pattern = "*.bam$", full.names = TRUE)
bam.files = bam.files[grep('tbx', bam.files, invert = TRUE)]

peak.list = list.files(path = DIR.peaks, pattern = "*.bed", full.names = TRUE)

# prepare peak regions and background
peak.list = peak.list[grep("pooled|random", peak.list)]

source("functions_chipSeq.R")
if(addBackground){
  peaks = merge.peaks.macs2(peak.list[grep('random', peak.list, invert = TRUE)], merge.dist = NULL)
  bgs = merge.peaks.macs2(peak.list[grep('random', peak.list)], merge.dist = NULL)
}

##########################################
# count reads within peaks and background
##########################################
source("functions_chipSeq.R")

if(run.featureCount.quantification){
  fc = quantify.signals.within.peaks(peaks, bam.list = bam.files, isPairedEnd = TRUE)
  fc.bgs = quantify.signals.within.peaks(bgs, bam.list = bam.files, isPairedEnd = TRUE)
  
  save(fc, fc.bgs, file = paste0(resDir, "counts_withinPeaksAndBackground_byfeatureCount.Rdata")) 
}else{
  load(file = paste0(resDir, "counts_withinPeaksAndBackground_byfeatureCount.Rdata"))
}

##########################################
# design matrix
##########################################
source("functions_chipSeq.R")
bams = fc$targets;
bams = gsub('X.Volumes.groups.cochella.jiwang.Projects.Ariane.tbx_transcriptional_priming.data.Bams.', '', bams)
#bams[grep('tbx_merged', bams)] = "../data/Bams/tbx_90min_merged.bam"
design.matrix = make.design.matrix.from.file.list(bams)

design.matrix$condition[which(design.matrix$condition=="60min")] = "90min"
design.matrix$condition[which(design.matrix$condition=="140min")] = "200min"
design.matrix$factor[grep("Aba", design.matrix$factor)] = "ABa"
design.matrix$factor[grep("Abp", design.matrix$factor)] = "ABp"
design.matrix$factor.condition = paste0(design.matrix$condition, "_", design.matrix$factor)
design.matrix$sampleID = sapply(design.matrix$file.name, function(x) unlist(strsplit(gsub('.bam', '', x), "_"))[3])
design.matrix$file.name = paste0(design.matrix$factor, "_", design.matrix$condition, "_", design.matrix$sampleID)
design.matrix$batch = NA
design.matrix$batch[grep('738', design.matrix$file.name)] = 2
design.matrix$batch[grep('713', design.matrix$file.name)] = 1

stat = fc$stat
kk = which(stat$Status=="Assigned" | stat$Status== "Unassigned_NoFeatures")
ss = apply(stat[kk, -1], 2, sum)
names(ss) = design.matrix$file.name

annot.bg = fc.bgs$annotation 
annot.bg$GeneID = paste0('bg_', annot.bg$GeneID)
annot = rbind(fc$annotation, annot.bg)

index.peaks = grep("bg_", annot$GeneID, invert = TRUE)

counts = rbind(fc$counts, fc.bgs$counts)
colnames(counts) = design.matrix$file.name

pdfname = paste0(resDir, "Data_Qulity_Assessment_atacPeakSignals", version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 12)

source("functions_chipSeq.R")
kk = c(which(colnames(design.matrix) == "factor.condition"), which(colnames(design.matrix)=="batch"))
#kk = which(colnames(design.matrix) == "factor.condition")
res = DB.analysis(counts[index.peaks, ], design.matrix[, kk], Threshold.read.counts = 100, batch = TRUE, size.factors = ss/median(ss))

dev.off()

##########################################
# calculate the rpkm and remove the batch difference
##########################################
dds = DESeqDataSetFromMatrix(counts[index.peaks,], DataFrame(design.matrix), design = ~ factor.condition)
#dds <- dds[ rowSums(counts(dds)) > 20, ]
sizeFactors(dds) = sizeFactors(res)
fpm = fpm(dds, robust = TRUE)
rownames(fpm) = annot$GeneID[index.peaks]

fpkm = fpm/annot$Length[index.peaks]*10^3

xx = log2(fpkm)
xx[which(rownames(xx)=='chrV_10647106_10647681'), ]

if(batch.removal){
  
  cat('remove the batch effect using limma \n')
  require('limma')
  #cat('remove the batch effect using ComBat \n')
  #require("sva")
  
  jj = which(design.matrix$factor == 'ABa' | design.matrix$factor == 'ABp')
  
  pdfname = paste0(resDir, "batchCorrection_limma", version.analysis, ".pdf")
  pdf(pdfname, width = 16, height = 12)
  
  source("functions_chipSeq.R")
  bc = remove.batch.atacseq(xx[,jj], design.matrix[jj, ], method = 'limma')
  
  dev.off()
  
}

# peak2gene assignment
if(run.peak2geneAssignment){
  peak.coord = data.frame(t(sapply(rownames(xx), function(x) unlist(strsplit(gsub("bg_", "",x), "_")))))
  
  source('functions_chipSeq.R')
  res.wsize.2kb = customizedAssignment.peak2gene(peak.coord = peak.coord, window.size=2000, annotation='wormbase')
  
  res.wsize.2kb = data.frame(peak.coord, res.wsize.2kb, stringsAsFactors = FALSE)
  rownames(res.wsize.2kb) = rownames(peak.coord)
  
  save(res.wsize.2kb, file = paste0(resDir, "peak2gene_assignment_4atacPeaks.Rdata")) 
}else{
  load(file = paste0(resDir, "peak2gene_assignment_4atacPeaks.Rdata")) 
}
 
# reorder the samples
o1 = order(design.matrix$factor, as.numeric(gsub("min", '',design.matrix$condition)))
time0 = grep('sorted.2to8cell', design.matrix$factor)
o1 = c(time0, setdiff(o1, time0) )
design.matrix = design.matrix[o1, ]
xx = xx[, o1]

#save(pp, res.wsize.2kb, file=paste0('Rdata/Merged_Peaks_macs2_p_5_filtered_N2', version.analysis, 'gene_assignment_WBcel235.Rdata'))
#res = res.wsize.2kb;
res = data.frame(res.wsize.2kb, xx, stringsAsFactors = FALSE)
colnames(res)[c(1:3)] = c('peak.chr', 'peak.start', 'peak.end') 

#kk = grep('lsy-6', res$gene);
res[which(rownames(res)=='chrV_10647106_10647681'),  grep('min_', colnames(res))]
save(res, design.matrix, file = paste0(resDir, "atacPeakSignals_geneAssignment_allReplicates.Rdata"))

write.csv(res, file = paste0(tableDir, "normalized_rpkm_for_atacSeqPeaks_background_geneAssignment.txt"),
            col.names = TRUE, row.names = TRUE, quote = FALSE)


########################################################
########################################################
# Section : clustering or pairwise comparisons using DESeq2 (not used here)
# due to the complex design (lineage + time points)
########################################################
########################################################
run.pairwise.Comparison.DESeq2 = FALSE
determine.groups.by.clustering = TRUE

if(run.pairwise.Comparison.DESeq2){
  dds = estimateDispersions(dds, fitType = "parametric")
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  plotDispEsts(dds, ylim=c(0.001, 10), cex=0.6)
  
  dds = nbinomWaldTest(dds, betaPrior = TRUE)
  resultsNames(dds)
  
  res1 <- results(dds, contrast = c("conds", "UNC3866", "Negative.Control.UNC4219"));
  res2 = results(dds, contrast = c("conds", "UNC4976", "Negative.Control.UNC4219"));
  summary(res1)
  res1 = as.data.frame(res1);
  summary(res2)
  res2 = as.data.frame(res2);
  
  write.table(fpm, file = paste0(tableDir, "normalized_readCounts_for_", prot, version.analysis, ".txt"), sep = "\t",
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(res1, file = paste0(tableDir, "DB_analysis_using_DESeq2_for_", prot, "_UNC3866_vs_Control",  version.analysis, ".txt"), sep = "\t",
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(res2, file = paste0(tableDir, "DB_analysis_using_DESeq2_for_", prot, "_UNC4976_vs_Control", version.analysis, ".txt"), sep = "\t",
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  
}else{
  
  load(file = paste0(resDir, "atacPeakSignals_geneAssignment_allReplicates.Rdata"))
  
  if(determine.groups.by.clustering){
    
    pks = as.matrix(res[, grep('min_', colnames(res))])
    kk = apply(pks, 1, function(x) all(x == -Inf) )
    pks = pks[!kk, ]
    
    pks[which(pks==-Inf)] = 0
    
    pks[which(rownames(pks)=='chrV_10647106_10647681'),  ]
    pks = pks[, -1]
    pks[which(rownames(pks)=='chrV_10647106_10647681'),  ]
    design.matrix = design.matrix[match(colnames(pks), design.matrix$file.name), ]
    
    pdfname = paste0(resDir, "clustering_atacseqPeaks_preprosessingWithLimma", version.analysis, ".pdf")
    pdf(pdfname, width = 16, height = 12)
    
    source('functions_chipSeq.R')
    #pks = merge.biologicalReplicates(pks)
    #pks[which(rownames(pks)=='chrV_10647106_10647681'),  ]
    #source('functions_chipSeq.R')
    cor.cutoff = 0.6 
    
    groups = clustering.peak.signals(pks, design.matrix, sd.cutoff = 0.7, cor.cutoff = 0.6, plot.grouping.result = TRUE)
    
    dev.off()
    
    res = data.frame(res[, c(1:11)],  pks, groups, stringsAsFactors = FALSE)
    
    load(file = paste0(resDir, "atacPeakSignals_geneAssignment_clustered.Rdata"))
    write.table(res, file = paste0(tableDir, "normalized_rpkm_for_atacSeqPeaks_background_geneAssignment_clustered.txt"),
                sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
    save(res, file = paste0(resDir, "atacPeakSignals_geneAssignment_clustered.Rdata"))
    
    
  }else{ # test a trendy package design for time series; but does not work well
    library(Trendy)
    require(stats)
    pks = as.matrix(res[, grep('min_', colnames(res))])
    #time.vector = 1:ncol(pks)
    time.vector <- c(1, rep(c(2:7), each=2))
    names(time.vector) = colnames(pks)
    
    xx <- trendy(Data = pks, tVectIn = time.vector, maxK=5, NCores = 6, minNumInSeg = 1, meanCut = 0)
    
    yy <- results(xx)
    yy.top <- topTrendy(yy, adjR2Cut = 0.3)
    yy.top$AdjustedR2
    #res.trend2 <- trendHeatmap(res.top2)
  }
  
}


########################################################
########################################################
# Section : small analysis for some details 
# 
########################################################
########################################################

##########################################
# ovlerapping between tbx peaks and the peak cluster where lsy-6 is
##########################################
load(file = paste0(resDir, "atacPeakSignals_geneAssignment_clustered.Rdata"))
cluster.sel = unique(res$groups[grep('lsy-6', res$gene)])

df = data.frame(res[which(res$groups == cluster.sel), c(1:3)], stringsAsFactors = FALSE)

xx = makeGRangesFromDataFrame(df)  # strand value "." is replaced with "*"

yy = readPeakFile(peak.files[grep('tbx_90min', peak.files)])
if(seqlevelsStyle(yy) != "UCSC") seqlevelsStyle(yy) = "UCSC";

length(xx)
length(xx[overlapsAny(xx, yy)])
length(countOverlaps(xx,yy))
sum(countOverlaps(xx,yy)>0)
