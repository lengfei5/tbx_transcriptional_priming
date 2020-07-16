##########################################################################
##########################################################################
# Project: tbx priming 
# Script purpose: atac-seq data analysis for paper revision
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Jul 15 12:38:38 2020
##########################################################################
##########################################################################
peak.Dir = c('../../all_ATACseq_for_Manuscript/Peaks/macs2')

peak.files = list.files(path = peak.Dir, pattern = "*.xls", full.names = TRUE)
#peak.files = c(peak.files, list.files(path = "../../R6729_atac/Peaks/macs2", 
#                        pattern = "*.xls", full.names = TRUE))
# import tbx peaks
# peak.files = c(peak.files, "../data/tbx_90min_peaks_merged_macs2_p_5_filtered_N2_gene_assignment_TSS_WBcel235_analysis_v2.bed")

outDir = paste0("../results/paper_revision")
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
#library("DiffBind")
#library("IRanges");
source('functions_chipSeq.R')

###############################
# make design matrix 
###############################
source('functions_chipSeq.R')
design.matrix = make.design.matrix.from.file.list(peak.files = peak.files, varnames = c('genetic.background', 'lineage', 'time.point'))

if(manual.modifySampleInfos){
  design.matrix$time.point[which(design.matrix$time.point=="60min")] = "90min"
  design.matrix$time.point[which(design.matrix$time.point=="140min")] = "200min"
  design.matrix$time.point[which(design.matrix$time.point=="220min")] = "200min"
  design.matrix$time.point[which(design.matrix$time.point=="350min")] = "330min"
  design.matrix$time.point[grep('batch|ASE', design.matrix$time.point)] = "sorted.ASE"
  # design.matrix$condition[which(design.matrix$condition=="140min")] = "200min"
  # design.matrix$factor[grep("Aba", design.matrix$factor)] = "ABa"
  # design.matrix$factor[grep("Abp", design.matrix$factor)] = "ABp"
  # design.matrix$factor.condition = paste0(design.matrix$condition, "_", design.matrix$factor)
  # design.matrix = design.matrix[order(design.matrix$condition, design.matrix$factor), ]
}

##########################################
# Import and merge peak files
##########################################
cat("-- import peak files as GRanges objects \n")
peaks.list = c()
peaks.all = NULL
design.matrix$nb.peaks = NA

for(n in 1:nrow(design.matrix)){
  #n = 14
  cat(n, '--', design.matrix$file.name[n], design.matrix$file.path[n], '\n')
  xx = readPeakFile(design.matrix$file.path[n], as = "GRanges");
  cat(length(xx), ' peaks found \n')
  if(seqlevelsStyle(xx) != "UCSC") seqlevelsStyle(xx) = "UCSC";
  
  design.matrix$nb.peaks[n] = length(xx)
  
  peaks.list = c(peaks.list, xx)
  if(n==1) {
    peaks.all = xx; 
  }else{
    peaks.all = union(peaks.all, xx, ignore.strand = TRUE)
  }
}
names(peaks.list) = design.matrix$file.name

lsy6.peaks = data.frame(chr = "chrV", start = c(10646957, 10647225), 
                        end = c(10647220, 10647697), score = c(1000, 1000))
lsy6.peaks = makeGRangesFromDataFrame(lsy6.peaks)

#peaks.without.lsy6.peaks = subsetByOverlaps(peaks.all, lsy6.peaks, invert = FALSE)
index.overlap = findOverlaps(peaks.all, lsy6.peaks)
peaks.without.lsy6.peaks = peaks.all[-index.overlap@from]
#subsetByOverlaps(lsy6.peaks, peaks.all)

peaks.xx = c(peaks.without.lsy6.peaks, lsy6.peaks)

peaks.all = reduce(peaks.xx)

subsetByOverlaps(peaks.all, lsy6.peaks, invert = FALSE)

save(design.matrix, peaks.all, peaks.list, file = paste0(outDir, '/design.matrix_peak.list_pooledpeaks.Rdata'))
export(peaks.all, con = paste0(outDir,  "/ABa_ABp_ASE_pooledPeaks.bed"))

########################################################
########################################################
# Section: quantify signals within peaks and background regions
# Differential Binding analysis
########################################################
########################################################
DIR.bams = "../../all_ATACseq_for_Manuscript/alignments/BAMs_unique_rmdup"
DIR.peaks = outDir

resDir = outDir
tableDir = paste0(resDir, "/tables/")
version.analysis = "_20200715"

addBackground = FALSE
run.featureCount.quantification = TRUE
run.peak2geneAssignment = FALSE

batch.removal = FALSE

if(!dir.exists(resDir)) system(paste0('mkdir -p ', resDir))
if(!dir.exists(tableDir)) system(paste0('mkdir -p ', tableDir))

bam.files = list.files(path = DIR.bams, pattern = "*.bam$", full.names = TRUE)
bam.files = bam.files[grep('tbx', bam.files, invert = TRUE)]

peak.list = list.files(path = DIR.peaks, pattern = "*.bed", full.names = TRUE)

# prepare peak regions and background
#peak.list = peak.list[grep("pooled|random", peak.list)]

source("functions_chipSeq.R")
if(addBackground){
  peaks = merge.peaks.macs2(peak.list[grep('random', peak.list, invert = TRUE)], merge.dist = NULL)
  bgs = merge.peaks.macs2(peak.list[grep('random', peak.list)], merge.dist = NULL)
}else{
  peaks = merge.peaks.macs2(peak.list, merge.dist = NULL)
}

# double check the lsy-6 peaks
kk = which(peaks$seqnames == 'chrV' & peaks$start >10646950 & peaks$end < 10647698)
peaks[kk, ]

##########################################
# count reads within peaks and background
##########################################
source("functions_chipSeq.R")

if(run.featureCount.quantification){
  fc = quantify.signals.within.peaks(peaks, bam.list = bam.files, isPairedEnd = TRUE)
  if(addBackground){
    fc.bgs = quantify.signals.within.peaks(bgs, bam.list = bam.files, isPairedEnd = TRUE)
    save(fc, fc.bgs, file = paste0(resDir, "/counts_withinPeaksAndBackground_byfeatureCount.Rdata")) 
  }else{
    save(fc, file = paste0(resDir, '/counts_withinPeaks_byfeatureCount.Rdata')) 
  }
}

##########################################
# process counts within peaks
##########################################
load(file = paste0(resDir, '/counts_withinPeaks_byfeatureCount.Rdata')) 

source("functions_chipSeq.R")

bams = basename(bam.files)
bb = gsub('_', '.', bams)
mm = match(bb, fc$targets)
bams = bams[mm]

##########################################
# !!!!! correct sample info due to sample swapping in R6729 
##########################################
kk1 = intersect(grep('MLC1480_ABa', bams), grep('738', bams))
kk2 = intersect(grep('MLC1480_ABp', bams), grep('738', bams))
bams[kk1] = sapply(bams[kk1], function(x) gsub('ABa', 'ABp', x))
bams[kk2] = sapply(bams[kk2], function(x) gsub('ABp', 'ABa', x))

design.matrix = make.design.matrix.from.file.list(peak.files = bams, 
                                                  varnames = c('genetic.background', 'lineage', 'time.point', 'sampleID'))

if(manual.modifySampleInfos){
  design.matrix$time.point[which(design.matrix$time.point=="60min")] = "90min"
  design.matrix$time.point[which(design.matrix$time.point=="140min")] = "200min"
  design.matrix$time.point[which(design.matrix$time.point=="220min")] = "200min"
  design.matrix$time.point[which(design.matrix$time.point=="350min")] = "330min"
  design.matrix$time.point[grep('batch|ASE', design.matrix$time.point)] = "sorted.ASE"
  
}

stat = fc$stat
kk = which(stat$Status=="Assigned" | stat$Status== "Unassigned_NoFeatures")
ss = apply(stat[kk, -1], 2, sum)
names(ss) = design.matrix$file.name

#annot.bg = fc.bgs$annotation 
#annot.bg$GeneID = paste0('bg_', annot.bg$GeneID)
annot = fc$annotation

counts = fc$counts
design.matrix$condition = paste0(design.matrix$genetic.background, '_', design.matrix$lineage, '_', design.matrix$time.point)
design.matrix$file.name = paste0(design.matrix$condition, '_', design.matrix$sampleID)

colnames(counts) = design.matrix$file.name
rownames(counts) = annot$GeneID
rownames(design.matrix) = design.matrix$file.name
#index.peaks = grep("bg_", annot$GeneID, invert = TRUE)

##########################################
# Normalization and comparison 
##########################################
source("functions_chipSeq.R")

dds = DESeqDataSetFromMatrix(counts, DataFrame(design.matrix), design = ~ condition)
#dds <- dds[ rowSums(counts(dds)) > Threshold.read.counts, ]
dds <- dds[ rowSums(counts(dds)) > 10, ]
sizeFactors(dds) = ss/median(ss)
#dds <- estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)


#fpkm = fpm/annot$Length[index.peaks]*10^3
xx = log2(fpm + 2^-6)
index.lsy6 = c(which(rownames(xx)=='chrV_10646957_10647220'), which(rownames(xx) == 'chrV_10647225_10647697'))
xx[index.lsy6, grep('90min*', colnames(xx))]
xx[index.lsy6, ]

##########################################
# merge biological replicates 
##########################################
cc = unique(design.matrix$condition)
yy = matrix(NA, ncol = length(cc), nrow = nrow(xx))
rownames(yy) = rownames(xx)
colnames(yy) = cc
for(n in 1:length(cc))
{
  cat(n, '--', cc[n], '\n')
  kk = which(design.matrix$condition == cc[n])
  print(colnames(xx)[kk])
  yy[,n] = apply(xx[,kk], 1, median)
}

##########################################
# plots for paper revision
##########################################
pdfname = paste0(resDir, "/atacPeakSignals_peakSignals_analysis", version.analysis, ".pdf")
pdf(pdfname, width = 12, height = 10)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

## pca for overview
pca = plotPCA(vsd, intgroup = 'condition', returnData=FALSE)
plot(pca)

pca.saved = plotPCA(vsd, intgroup = 'condition', returnData=TRUE)
mm = match(rownames(pca.saved), design.matrix$file.name)
pca.saved = data.frame(pca.saved, design.matrix[mm, ])
pca.saved$lineage.genotype = paste0(pca.saved$lineage, '_', pca.saved$genetic.background)
ggp = ggplot(data=pca.saved, aes(PC1, PC2, shape=time.point, color = lineage.genotype )) + 
  geom_point(size=3) 
#+
  #scale_shape_manual(values=c(21:24)) + 
  #scale_alpha_manual(values=c("MLC1480"=0, " MLC2309"=1, 'MLC2310' = 2, 'otls252.253' = 3)) 
  #+
  #geom_text(hjust = 0.7, nudge_y = 2.5, size=2.5)
plot(ggp);

## sample correlations
library(corrplot)
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                           "cyan", "#007FFF", "blue","#00007F"))

sels = c("MLC1480_ABa_90min", "MLC1480_ABp_90min", "MLC1480_ABa_200min", "MLC1480_ABp_200min",
         "MLC2309_ABa_90min", "MLC2309_ABp_90min", "MLC2309_ABa_200min", "MLC2309_ABp_200min")
M <- cor(yy[, match(sels, colnames(yy))], use = 'everything')
#corrplot(M, method="circle", type = 'upper', order="hclust")
corrplot(M, method="number", type = 'upper')

M = cor(yy[, grep('ASE', colnames(yy))], use = 'everything')
corrplot(M, method="number", type = 'upper')
#corrplot(M, method="ellipse", order="hclust", tl.cex=1.2, cl.cex=0.7, tl.col="black", 
#         addrect=ceiling(ncol(yy)/2), col=col1(100), rect.col=c('green'), rect.lwd=2.0)


## scatter plots 
compares = list(c("MLC1480_ABa_90min", "MLC1480_ABp_90min"),
                c("MLC1480_ABa_200min", "MLC1480_ABp_200min"), 
                c("otls252.253_ASEL_sorted.ASE", "otls252.253_ASER_sorted.ASE"),
                c("MLC2309_ABa_90min", "MLC1480_ABa_90min"),
                c("MLC2309_ABp_90min", "MLC1480_ABp_90min"),
                c("MLC2310_ASE_sorted.ASE", "otls252.253_ASEL_sorted.ASE"), 
                c("MLC2310_ASE_sorted.ASE", "otls252.253_ASER_sorted.ASE")
                )

for(n in 1:length(compares))
{
  #n = 1
  kk = match(compares[[n]], colnames(yy))
  plot(yy[,kk], cex= 0.5)
  abline(0, 1, lwd=2.0, col='red')
  #points(t(yy[index.lsy6[1], kk]), col='orange', cex= 1.5, pch=16)
  points(t(yy[index.lsy6[2], kk]), col='blue', cex= 2.0, pch=16)
  #points(yy[index.lsy6, kk], col='blue', cex= 1.5, pch=16)
  
}

dev.off()

# ##########################################
# # calculate the rpkm and remove the batch difference
# ##########################################
# dds = DESeqDataSetFromMatrix(counts[index.peaks,], DataFrame(design.matrix), design = ~ factor.condition)
# #dds <- dds[ rowSums(counts(dds)) > 20, ]
# sizeFactors(dds) = sizeFactors(res)
# fpm = fpm(dds, robust = TRUE)
# rownames(fpm) = annot$GeneID[index.peaks]
# 
# fpkm = fpm/annot$Length[index.peaks]*10^3
# 
# xx = log2(fpkm)
# xx[which(rownames(xx)=='chrV_10647106_10647681'), ]
# 
# # 
# # if(batch.removal){
# #   
# #   cat('remove the batch effect using limma \n')
# #   require('limma')
# #   #cat('remove the batch effect using ComBat \n')
# #   #require("sva")
# #   
# #   jj = which(design.matrix$factor == 'ABa' | design.matrix$factor == 'ABp')
# #   
# #   pdfname = paste0(resDir, "batchCorrection_limma", version.analysis, ".pdf")
# #   pdf(pdfname, width = 16, height = 12)
# #   
# #   source("functions_chipSeq.R")
# #   bc = remove.batch.atacseq(xx[,jj], design.matrix[jj, ], method = 'limma')
# #   
# #   dev.off()
# #   
# # }


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
