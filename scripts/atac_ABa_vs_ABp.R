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
run.featureCount.quantification = TRUE
run.pairwise.Comparison.DESeq2 = FALSE
batch.removal = FALSE


if(!dir.exists(resDir)) system(paste0('mkdir -p ', resDir))
if(!dir.exists(tableDir)) system(paste0('mkdir -p ', tableDir))

bam.files = list.files(path = DIR.bams, pattern = "*.bam$", full.names = TRUE)
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
  fc = quantify.signals.within.peaks(peaks, bam.list = bam.files)
  fc.bgs = quantify.signals.within.peaks(bgs, bam.list = bam.files)
  save(fc, fc.bgs, file = paste0(resDir, "counts_withinPeaksAndBackground_byfeatureCount.Rdata")) 
}else{
  load(file = paste0(resDir, "counts_withinPeaksAndBackground_byfeatureCount.Rdata"))
}

bams = bam.files;
bams[grep('tbx_merged', bams)] = "../data/Bams/tbx_90min_merged.bam"  
design.matrix = make.design.matrix.from.file.list(bams)
design.matrix$condition[which(design.matrix$condition=="60min")] = "90min"
design.matrix$condition[which(design.matrix$condition=="140min")] = "200min"
design.matrix$factor[grep("Aba", design.matrix$factor)] = "ABa"
design.matrix$factor[grep("Abp", design.matrix$factor)] = "ABp"
design.matrix$factor.condition = paste0(design.matrix$condition, "_", design.matrix$factor)
design.matrix$sampleID = sapply(design.matrix$file.name, function(x) unlist(strsplit(gsub('.bam', '', x), "_"))[3])
design.matrix$file.name = paste0(design.matrix$factor, "_", design.matrix$condition, "_", design.matrix$sampleID)

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
kk = which(colnames(design.matrix) == "factor.condition")
res = DB.analysis(counts[index.peaks, ], design.matrix[, kk], Threshold.read.counts = 20)

dev.off()

##########################################
# calculate the rpkm and remove the batch difference
##########################################
dds = DESeqDataSetFromMatrix(counts, DataFrame(design.matrix), design = ~ factor.condition)
#dds <- dds[ rowSums(counts(dds)) > 20, ]
sizeFactors(dds) = sizeFactors(res)
fpm = fpm(dds, robust = TRUE)
rownames(fpm) = annot$GeneID

fpkm = fpm/annot$Length*10^3

xx = log2(fpkm)

# reorder the samples
o1 = order(as.numeric(gsub("min", '',design.matrix$condition)), design.matrix$factor)
design.matrix = design.matrix[o1, ]
xx = xx[, o1]


# peak2gene assignment
makeGRangesFromDataFrame(df)
load(file= paste0('Rdata/Merged_Peaks_macs2_p_5_filtered_N2', version.analysis,'.Rdata'))
#load(file='Rdata/Merged_Peaks_macs2_p_5_filtered_N2.Rdata')
pp = mergedpeaks;

source('functions_chipSeq.R')
#res.wsize.1kb = PeakAnnotation.refseq(pp, 1000)
res.wsize.2kb = PeakAnnotation.customized(pp, window.size=2000, annotation='wormbase')


#save(pp, res.wsize.2kb, file=paste0('Rdata/Merged_Peaks_macs2_p_5_filtered_N2', version.analysis, 'gene_assignment_WBcel235.Rdata'))

#load(file='Rdata/Merged_Peaks_macs2_p_5_filtered_N2_gene_assignment_WBcel235.Rdata')
names =  c('chr.peak', 'start.peak', 'end.peak', 'width.peak', 'strand.peak', 
           'chr.gene', 'start.gene', 'end.gene', 'width.gene', 'strand.gene', 
           'WormBase.Gene.ID', 'gene',  'window.size')
#colnames(res.wsize.1kb) = names
colnames(res.wsize.2kb) = names

res = res.wsize.2kb;
res = data.frame(res, stringsAsFactors = FALSE)
res$chr.peak = as.character(res$chr.peak)
MakePeakNames = function(x){paste(x, sep='', collapse = '_')}
peaknames = rep(NA, nrow(res))
for(n in 1:nrow(res))
{
  #n = 2
  peaknames[n] = MakePeakNames(res[n, c(1:3)])
}

res = data.frame(peaknames, res, stringsAsFactors = FALSE)
colnames(res)[1] = 'peak.name'

Test.peaks.association = FALSE
if(Test.peaks.association)
{
  kk = grep('lsy-6', res$gene);
  mm = match(res$peak.name, unique(res$peak.name[kk]))
  res[which(!is.na(mm)==TRUE), ]
}

if(batch.removal){
  cat('remove the batch effect using ComBat \n')
  require("sva")
  
  design.matrix$batch = NA
  design.matrix$batch[grep('738', design.matrix$file.name)] = 2
  design.matrix$batch[grep('713', design.matrix$file.name)] = 1
  
  jj = which(design.matrix$batch==1 | design.matrix$batch == 2)
  batch = design.matrix$batch[jj];
  design.tokeep = design.matrix[jj,]
  mod = model.matrix(~ as.factor(factor.condition), data = design.tokeep);
  yy = ComBat(dat=xx[,jj], batch=batch, mod=mod, par.prior=TRUE, ref.batch = 2)
}



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
  
  #pdfname = paste0(resDir, "results_DB_analysis_ChIPseq_", prot, version.analysis, ".pdf")
  #pdf(pdfname, width = 12, height = 10)
  
  # plot(res1$, res2$log2FoldChange)
  #plot(apply(fpm[, c(1, 2)], 1, mean), apply(fpm[, c(3, 4)], 1, mean), log='xy', cex=0.7, xlab='control', ylab= "UNC3866");
  #abline(0, 1, col='red', lwd=2.0)
  #plot(apply(fpm[, c(1, 2)], 1, mean), apply(fpm[, c(5, 6)], 1, mean), log='xy', cex=0.7, xlab= "control", ylab="UNC4976")
  #abline(0, 1, col='red', lwd=2.0)
  
  #dev.off();
  
  write.table(fpm, file = paste0(tableDir, "normalized_readCounts_for_", prot, version.analysis, ".txt"), sep = "\t",
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(res1, file = paste0(tableDir, "DB_analysis_using_DESeq2_for_", prot, "_UNC3866_vs_Control",  version.analysis, ".txt"), sep = "\t",
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(res2, file = paste0(tableDir, "DB_analysis_using_DESeq2_for_", prot, "_UNC4976_vs_Control", version.analysis, ".txt"), sep = "\t",
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  
}


