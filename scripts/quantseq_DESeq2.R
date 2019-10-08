##########################################################################
##########################################################################
# Project:
# Script purpose: run DESeq2 normalization on quant-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Apr  2 11:41:17 2019
##########################################################################
##########################################################################
Rfuncitons = '../../../../scripts/functions/RNAseq_functions.R'

version.Data = 'quantseq_ABa_ABp_R7138';
version.analysis = paste0("_", version.Data, "_20180402")

### Directories to save results 
#design.file = "../exp_design/Libaries_time_series_spikIns.xlsx"
dataDir = "../data/quantseq/featurecounts_R7183"

resDir = "../results/quantseq"
tabDir =  paste0(resDir, "/tables/")
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

##################################################
##################################################
## Section: import design matrix and prepare count table
##################################################
##################################################
Processing.design.matrix = TRUE
if(Processing.design.matrix){
  design = data.frame(c('80391', '80386', '80387', '80388', '80392'), 
                      c('16', '60', '60', '140', '140'), c('ABa.ABp', 'ABp', 'ABa', 'ABp', 'ABa'), stringsAsFactors = FALSE)
  colnames(design) = c('SampleID', 'time', 'lineage')
  design = design[, c(1, 3, 2)]
}

## make the data table
xlist<-list.files(path=paste0(dataDir), pattern = "*_gene.featureCounts.txt$", full.names = TRUE) ## list of data set to merge

source(Rfuncitons)
all = cat.countTable(xlist, countsfrom = 'featureCounts')

source(Rfuncitons)
xx = process.countTable(all=all, design = design, ensToGeneSymbol = TRUE)
all = xx;

save(design, all, file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))

####################
## QC for cpm 
####################
QC.for.cpm = TRUE
if(QC.for.cpm){
  load(file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))
  design$conds = paste0(design$lineage, "_", design$time)
  
  index.qc = c(4)
  
  source(Rfuncitons)
  
  pdfname = paste0(resDir, "/Data_qulity_assessment", version.analysis, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  Check.RNAseq.Quality(read.count=all[, -1], design.matrix = rep(1, nrow(design)))
  
  dev.off()
  
}

####################
## DESeq2 normalization
####################
require('DESeq2')
dds <- DESeqDataSetFromMatrix(all[, -1], DataFrame(design), design = ~ conds)
sels = rowSums(counts(dds)) > 10
ggs = all$gene[sels]
dds <- dds[sels, ]
res = fpm(dds, robust = TRUE)

#colnames(cpm) = paste0(colnames(cpm), ".cpm")
colnames(res) = paste0(colnames(res), ".DESeq2Norm")

res = data.frame(ggs, res, stringsAsFactors = FALSE)
colnames(res)[1] = 'gene'

write.csv(res, file = paste0(tabDir, paste0("DESeq2Normalized_signals", version.analysis, ".csv")), 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = TRUE)
