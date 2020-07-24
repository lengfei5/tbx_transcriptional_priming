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
# make design matrix for peak files
###############################
source('functions_chipSeq.R')
design.matrix = make.design.matrix.from.file.list(peak.files = peak.files, varnames = c('genetic.background', 'lineage', 'time.point'))

if(manual.modifySampleInfos){
  #design.matrix$time.point[which(design.matrix$time.point=="60min")] = "90min"
  #design.matrix$time.point[which(design.matrix$time.point=="140min")] = "200min"
  #design.matrix$time.point[which(design.matrix$time.point=="220min")] = "200min"
  #design.matrix$time.point[which(design.matrix$time.point=="350min")] = "330min"
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
# Section : define lineage specific peaks using called peaks
# define ABa- and ABp- specific peaks by pooling replicates and two time points (90 and 200 min)
########################################################
########################################################
library(rtracklayer)

design.matrix$condition = paste0(design.matrix$genetic.background, '_', design.matrix$lineage, 
                                 '_', design.matrix$time.point)

#sels = which(design.matrix$condition != "330min" & design.matrix$factor != 'tbx' & design.matrix$condition != "200min")
#sels = which(design.matrix$factor.condition == '90min_ABa' | design.matrix$factor.condition == '90min_ABp')
#design.sels = design.matrix[sels,  ]
#peaks.sels = peaks.list[sels]

conds = unique(design.matrix$condition)
peaks = c()
peaknames = conds

pval = 10
pool.peaks.for.biologicalRep = FALSE
Check.replicate.peakOverlap = FALSE

for(cc in conds){
  # fac = 'ABa'
  kk = which(design.matrix$condition == cc)
  peaks.merged = c()
  
  for(k in kk) 
  {
    # k = 1
    xx = readPeakFile(design.matrix$file.path[k], as = "GRanges");
    if(seqlevelsStyle(xx) != "UCSC") seqlevelsStyle(xx) = "UCSC";
    with.p.values = "X.log10.pvalue." %in% colnames(mcols(xx))
    p10 <- xx[mcols(xx)[,"X.log10.pvalue."] > pval]
    
    if(length(peaks.merged) == 0) {
      peaks.merged = reduce(p10, ignore.strand = TRUE)
    }else{
      if(pool.peaks.for.biologicalRep){
        peaks.merged = reduce(union(peaks.merged, p10, ignore.strand = TRUE), ignore.strand = TRUE)
      }else{
        peaks.merged = reduce(intersect(peaks.merged, p10, ignore.strand = TRUE), ignore.strand = TRUE)
      }
    }
    cat(cc, ' -- all peaks :', length(xx), ' -- peaks <10^-', pval, ' : ',  length(p10), ' -- merged peaks', length(peaks.merged), '\n')
  }
  peaks = c(peaks, peaks.merged)
  
}

names(peaks) = conds

if(Check.replicate.peakOverlap){
  pdf(paste0(outDir, "/Comparison_ATAC_peaks_for_ABa_ABp_replicateOverlapPeaks.pdf"),
      width = 12, height = 8)
  source("functions_chipSeq.R") 
  
  ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=c('ABa', 'ABp'), connectedPeaks="keepAll", main='ABa.vs.ABp')
  
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
  dev.off()
}

##########################################
# define peak group of interest 
##########################################
Define.Groups.peaks = FALSE
if(Define.Groups.peaks){
  
  system(paste0('mkdir -p ', paste0(outDir, '/peak_groups')))
  
  ## wt_ABa_90min peaks
  p1 = peaks[[which(names(peaks) == "MLC1480_ABa_90min")]]
  #export(p1, con = paste0(outDir, '/peak_groups/peaks_wt_ABa_90min.bed'))
  
  ## wt_ABa_90min peaks not overlapped by wt_ABp_90min
  p2 = peaks[[which(names(peaks) == "MLC1480_ABp_90min")]]
  p11 = p1[!overlapsAny(p1, p2)]
  export(p11, con = paste0(outDir, '/peak_groups/peaks_wt_ABa_90min_not_sharedby_ABp.bed'))
  
  ## peaks for ASER and ASEL
  p3 = peaks[[which(names(peaks) == "otls252.253_ASEL_sorted.ASE")]]
  p4 = peaks[[which(names(peaks) == "otls252.253_ASER_sorted.ASE")]]
  
  p33 = p3[!overlapsAny(p3, peaks[[which(names(peaks) == "MLC1480_ABa_330min")]])]
  p44 = p4[!overlapsAny(p4, peaks[[which(names(peaks) == "MLC1480_ABp_330min")]])]
  
  export(p3, con = paste0(outDir, '/peak_groups/peaks_wt_ASEL.bed'))
  export(p4, con = paste0(outDir, '/peak_groups/peaks_wt_ASER.bed'))
  
  export(p33, con = paste0(outDir, '/peak_groups/peaks_wt_ASEL_vs_ABa.330min.bed'))
  export(p44, con = paste0(outDir, '/peak_groups/peaks_wt_ASER_vs_ABp.330min.bed'))
  
  
  ## group of genes similar to lsy-6: tbx-binding, open in ABa_wt_90min but closed in ABp_wt_90min
  p1 = peaks[[which(names(peaks) == "MLC1480_ABa_90min")]]
  p2 = peaks[[which(names(peaks) == "MLC1480_ABp_90min")]]
  p11 = p1[!overlapsAny(p1, p2)]
  
  tbx = "../data/tbx_90min_peaks_merged_macs2_p_5_filtered_N2_gene_assignment_TSS_WBcel235_analysis_v2.bed"
  tbx = readPeakFile(tbx, as = "GRanges");
  if(seqlevelsStyle(tbx) != "UCSC") seqlevelsStyle(tbx) = "UCSC";
  
  p11 = p11[overlapsAny(p11, tbx)]
  
  saveRDS(p11, file = paste0(outDir, '/ATACseq_peaks_similar_to_lsy6.rds'))
  
  subsetByOverlaps(p11, lsy6.peaks)
  export(p11, con = paste0(outDir, 'ATACseq_peaks_similar_to_lsy6.bed'))
  #lsy6.peaks
  
  # p.all = union(peaks[[1]], peaks[[2]], ignore.strand = TRUE)
  # p.all = union(p.all, peaks[[3]], ignore.strand = TRUE)
  # p.all = reduce(p.all, ignore.strand = TRUE)
  # 
  # export(p.all, con = paste0("../results/peakGroups/early_ABa_ABp_pooledPeaks.bed"))
  # export(peaks[[1]], con = paste0("../results/peakGroups/early_allPeaks.bed"))
  # export(peaks[[2]], con = paste0("../results/peakGroups/ABa_allPeaks.bed"))
  # export(peaks[[3]], con = paste0("../results/peakGroups/ABp_allPeaks.bed"))
  
}

########################################################
########################################################
# Section: quantify signals within peaks and background regions
# Differential Binding analysis
########################################################
########################################################
DIR.bams = "../../all_ATACseq_for_Manuscript/alignments/BAMs_unique_rmdup_old"
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

##########################################
# prepare peak regions and background
##########################################
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

##########################################
# scatter plots with highlighted peaks
##########################################
library(Signac)

peak.sels = c('chrV_10646957_10647220', 
              'chrV_10647225_10647697',
              'chrII_10316291_10316507',
              'chrV_20825034_20825321',
              'chrV_8387150_8387366',
              'chrV_14867871_14868107')
peak.sels = StringToGRanges(peak.sels, sep = c('_', '_'))

regions = StringToGRanges(rownames(yy), sep = c("_", "_"))
regions.sel = regions[overlapsAny(regions, peak.sels, minoverlap = 20L)]

gene.sels = c('gcy-5','gcy-7', 'lsy-6.promoter', 'lsy-6.enhancer',  'gcy-14',
              'gcy-22')
index.sels = match(regions.sel, regions)

#index.lsy6 = c(which(rownames(xx)=='chrV_10646957_10647220'), which(rownames(xx) == 'chrV_10647225_10647697'))
#xx[index.lsy6, grep('90min*', colnames(xx))]
#xx[index.lsy6, ]

## scatter plots 
compares = list(c("MLC1480_ABa_90min", "MLC1480_ABp_90min"),
                c("MLC1480_ABa_200min", "MLC1480_ABp_200min"), 
                c("MLC1480_ABa_90min", "MLC2309_ABa_90min"),
                c("otls252.253_ASEL_sorted.ASE", "MLC2310_ASE_sorted.ASE"),
                c("otls252.253_ASEL_sorted.ASE", "otls252.253_ASER_sorted.ASE")
                
                )

pdfname = paste0(resDir, "/atacPeakSignals_peakSignals_comparison_fianl_v2", version.analysis, ".pdf")
pdf(pdfname, width = 12, height = 12)
par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(6,16,2,0.8)+0.1, tcl = -0.3)
par(pty="s")

for(n in 1:length(compares))
{
  # n = 3
  kk = match(compares[[n]], colnames(yy))
  plot(yy[,kk], cex= 0.2, xlim = c(-5, 15), ylim = c(-5, 15))
  abline(0, 1, lwd=2.0, col='red')
  points(yy[index.sels[3:4], kk], col='blue', cex= 1.2, pch=16)
  text(yy[index.sels[3:4], kk], gene.sels[3:4], offset = 0.5, pos = 4, cex = 0.8)
  rr = yy[,kk[1]] - yy[,kk[2]]
  if(n>=5){
    pch = c(15:18)
    points(yy[index.sels[c(1:2, 5:6)], kk], col='orange', cex= 1.2, pch=pch)
    #text(yy[index.sels[c(1:2, 5:6)], kk], gene.sels[c(1:2, 5:6)], offset = 1.0, pos = 4)
    legend('topleft', legend =  gene.sels[c(1:2, 5:6)], pch = pch, bty = 'n', col = 'orange')
  }
  
  #points(yy[index.lsy6, kk], col='blue', cex= 1.5, pch=16)
  
}

dev.off()


##########################################
# heatmap for genes or peaks similar to lsy-6 
##########################################
library(Signac)
pp = readRDS(file = paste0(outDir, '/ATACseq_peaks_similar_to_lsy6.rds'))

regions = StringToGRanges(rownames(yy), sep = c("_", "_"))
regions.sel = regions[overlapsAny(regions, pp, minoverlap = 100L)]

kk = match(regions.sel, regions)
mm = match(c("MLC1480_ABa_90min", "MLC1480_ABa_200min", "MLC1480_ABa_330min", 
             "MLC1480_ABp_90min", "MLC1480_ABp_200min", "MLC1480_ABp_330min"), colnames(yy))

keep = yy[kk, mm]

library("pheatmap");
library("RColorBrewer");

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(keep, 
         scale = 'row',
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         gaps_col = 3,
         fontsize_row = 6,
         clustering_distance_rows = 'correlation'
         #clustering_distance_cols = sampleDists,
         #col = colors
         )

require(ChIPseeker)
require(TxDb.Celegans.UCSC.ce11.ensGene)
txdb <- TxDb.Celegans.UCSC.ce11.ensGene
peakAnnoList <- annotatePeak(regions.sel, TxDb=txdb, tssRegion=c(-2000, 2000), verbose=FALSE)
peakannot = data.frame(peakAnnoList, stringsAsFactors = FALSE)

load(file = '/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.Rdata')

peakannot$genes = as.character(annot$Gene.name[match(peakannot$geneId, annot$Gene.stable.ID)])
kk = which(is.na(peakannot$genes))

peakannot$genes[kk] = as.character(peakannot$geneId[kk])

peakannot$genes[which(peakannot$seqnames=='chrV' & peakannot$end == '10647697')]

rownames(keep) = make.unique(peakannot$genes)
pheatmap(keep, 
         scale = 'row',
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         gaps_col = 3,
         fontsize_row = 5,
         clustering_distance_rows = 'correlation'
         #clustering_distance_cols = sampleDists,
         #col = colors
)
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
