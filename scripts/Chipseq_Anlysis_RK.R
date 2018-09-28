##################################################
##################################################
## Project: Jorge's project
## Script purpose: 
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Tue Mar 13 15:19:17 2018
##################################################
##################################################

##################################################
##################################################
## Section: library, data and process design matrix
##################################################
##################################################
library("ChIPseeker");
library("rtracklayer")
#library(UpSetR);library("ChIPpeakAnno")
library("Vennerable") #install.packages("Vennerable", repos="http://R-Forge.R-project.org", type="source")
library("ggplot2");
library("GenomicFeatures");
library("GenomicRanges")
#library("DiffBind");library("IRanges");
source("functions_chipSeq.R") ## these functions are from Thomas Burkard 

### Import peaks and design matrix
version.analysis = 'peaks_ZFP462_vs_HP1_Dnmt'
peakDir = "Peaks/macs2_broad"

#design.file = ""
outDir = "result_compare_ZFP462_otherChIPs/"
if(!dir.exists(outDir)) dir.create(outDir)

xlist<-list.files(path=peakDir, pattern = "*macs2_broad_fdr_0.1_peaks.xls", full.names = TRUE)
#xlist = xlist[grep('_H', xlist)]

peaks.list = c() 
for(n in 1:length(xlist)){
  cat(n, '\n')
  xx = readPeakFile(xlist[n], as = "GRanges"); 
  peaks.list = c(peaks.list, xx)
  #eval(parse(text = paste0("pp.", n, " = xx")))
}

names(peaks.list) = basename(xlist)

kk = grep('ZFP462', names(peaks.list))
bg0 = grep("ATAC|atac", names(peaks.list))

for(k in kk)
{
  xx = peaks.list[[k]];
  cat('before filter', length(xx), ' peaks\n');
  for(bg in bg0)
  {
    xx = xx[!overlapsAny(xx, peaks.list[[bg]])];
  }
  cat('after filter', length(xx), ' peaks\n');
  peaks.list[[k]] = xx
}

names(peaks.list) = sapply(names(peaks.list), function(x) paste0(unlist(strsplit(as.character(x), "_"))[c(1:2)], collapse = "_"))

## edit the pair-wise comparisons
#cc = list()

#cc = cbind(rep(20, 19), c(1:19))
#cc = rbind(cc, c(5, 11), c(17,15), c(17, 16), 
#           c(10, 11), c(10, 12), c(10, 13), c(10,14), c(10, 7))
pdf(paste0(outDir, "/Comparison_Peaks_overlapping_for_ZFP462_macs2_", version.analysis, ".pdf"), width = 12, height = 8)

for(n in 1:19)
{
  # n = 1
  kk = c(n, 20)
  peaks = c()
  #peaks10 = c()
  peaknames = names(peaks.list)[kk]
  #peaknames = unlist(sapply(peaknames, function(x) gsub("_macs2_broad_fdr_0.1_peaks.xls", "", x)))
  
  for(k in kk) 
  {
    p = peaks.list[[k]]
    #p = readPeakFile(ff$file.path[k], as = "GRanges");
    #eval(parse(text = paste0("p = pp.", k)));
    #p10 <- p[mcols(p)[,"X.log10.pvalue."] > pval];
    #p <- mergeWindows(p, tol=2000L, ignore.strand = TRUE)
    p = reduce(p); 
    #p10 = reduce(p10);
    peaks= c(peaks, p);
    #peaks10= c(peaks10, p10);
    #peaknames = c(peaknames, ff$file.name[k])
    #if(DB.Analysis){
    #  peaknames = c(peaknames, ff$sample.condition[k])
    #}else{
    #}
    #test = basename(macs2.files[k])
    #test = find.sample.names(test, n2=2)
    #test = unlist(strsplit(as.character(test), '_'))
    #test =test[length(test)]
    #nb.reads = bamstatistics[grep(test, bamstatistics[,2]), 5]
    #peaknames = c(peaknames, paste(ff$file.name[k], test,  sep='_', collapse = '_'))
  }
  
  ol.peaks <- makeVennDiagram(peaks, NameOfPeaks=peaknames, connectedPeaks="keepAll")
  v <- venn_cnt2venn(ol.peaks$vennCounts)
  try(plot(v))
  
}

dev.off()

pp = union(peaks.list[[19]], peaks.list[[20]], ignore.strand=TRUE)
xx = as.data.frame(pp)
df <- data.frame(seqnames=seqnames(pp), starts=start(pp), ends=end(pp), names=c(rep(".", length(pp))),
                                  scores=c(rep(".", length(pp))), strands=strand(pp))

sel.chr = paste0('chr', c(1:19, 'X', 'Y'));
mm = match(df$seqnames, sel.chr)
df = df[which(!is.na(mm)==TRUE), ]

write.table(df, file = paste0(outDir, "ZFP462_peaks_excluded_ATAC_mergedReplicates.bed"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

### select samples to check using design matrix
Import.Design.Matrix = FALSE
if(Import.Design.Matrix){
  ## improt design matrix
  design = read.delim(paste0("ChIPseq_design_matrix_all.txt"), header = TRUE, as.is = 2)
  index = c()
  for(n in 1:nrow(design))
  {
    jj = grep(design$SampleID[n], xlist)
    if(length(jj)==1){index = c(index, jj)
    }else{cat("NOT FOUND sample for ", design$sampleID[n], "\n")}
  }
  macs2.files = xlist[index]
  bname = basename(macs2.files) 
  design.matrix = data.frame(macs2.files, bname, design, stringsAsFactors = FALSE)
  colnames(design.matrix) = c('file.path', 'file.name', 'sampleID', 'condition', 'factor')
  #ff = data.frame(ff, sapply(bname, find.samples.conditions, ID='samples'))
  o1 = with(design.matrix, order(factor, condition))
  design.matrix = design.matrix[with(design.matrix, order(factor, condition)), ];
}else{
  ## make design matrix based on the peak files
  macs2.files = xlist;
  bname = basename(macs2.files)
  xx = c()
  yy = c()
  for(n in 1:length(bname))
  {
    #n = 1;
    test = bname[n];
    test = unlist(strsplit(as.character(test), "_"))
    xx = c(xx, test[1])
    #test = unlist(strsplit(as.character(test), "[.]"))[-1]
    #test = paste0(test, collapse = ".")
    yy = c(yy, test[2])
    #bb = sapply(bname, function(x) gsub("..", ".", x))
  }
  design = data.frame(yy, xx)
  
  #design = data.frame(sapply(bname, find.samples.conditions, ID='conditions'), sapply(bname, find.samples.conditions, ID='samples'), 
  #                    stringsAsFactors = FALSE)
  design.matrix = data.frame(macs2.files, bname, design, stringsAsFactors = FALSE)
  colnames(design.matrix) = c('file.path', 'file.name', 'condition', 'factor')
  design.matrix = design.matrix[with(design.matrix, order(factor, condition)), ];
  #condition= sapply(bname, find.sample.names)
  #macs2.files <- paste("/Volumes/groups/bell/jiwang/Projects/Jorge/ChIP_seq_Jorge/peakcalling/macs2/", xlist, sep='')
}

factor.condition = paste0(design.matrix$factor, '_', design.matrix$condition)
design.matrix = data.frame(design.matrix, factor.condition, stringsAsFactors = FALSE)

##################################################
##################################################
## Section: Quality control after peak calling by macs2
## check the peak overlapping 
##################################################
##################################################
## import macs peaks as Genomic Range object

pdf(paste0(outDir, "/Comparison_Peaks_All_Replicates_overlapping_macs2_", version.analysis, ".pdf"), width = 12, height = 8)

source("functions_chipSeq.R") #
#par(mfrow=c(1,2))
Comparison.overlapping.peaks(design.matrix, peaks.list, toCompare="factor.condition") 

dev.off()

##################################################
##################################################
## Section: Peak extension, Peak-to-gene assignment, save tables and annotation plots
##################################################
##################################################
ChIPseq.peak.gene.assignment = TRUE
Save.Peak.Annotation = TRUE;

Chr2save = paste0('chr', c(1:19, 'X', 'Y'));
chroms = Chr2save;
Merge.close.peaks = TRUE
merge.dist = 2000;

if(ChIPseq.peak.gene.assignment)
{
  library(csaw)
  library(AnnotationDbi)
  library(AnnotationHub)
  
  ## import gene annotation
  #txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene")
  #saveDb(txdb, file="refGene_mm10.sqlite")
  txdb <- loadDb("../refGene_mm10.sqlite")
  
  ah = AnnotationHub()
  #ce = ah[["AH52238"]]
  #mm <- ah[ah$species=="Mus musculus" & ah$rdataclass=='OrgDb']
  mm = ah[["AH57974"]]
  entrez2symbol <- unique(select(mm, keys=keys(mm, keytype="ENTREZID"), columns=c("SYMBOL"), keytype="ENTREZID"))
  
  #blacklist = readPeakFile(, as = "GRanges")
  pdfname = paste0(outDir, "/Peaks_Annotation_",version.analysis, ".pdf")
  pdf(pdfname, width = 10, height = 6)
  
  for(n in 1:nrow(design.matrix))
  {
    # n = 1 
    cat(design.matrix$file.name[n], "\n")
    
    #grep('H3K27me3', macs2.files)
    #n = 6;
    p = readPeakFile(design.matrix$file.path[n], as = "GRanges"); 
    #p10 <- p[mcols(p)[,"X.log10.pvalue."] > pval];
    cat(length(p), '\n')
    if(design.matrix$factor[n] == "H2AK119" | design.matrix$factor[n] == "H2AUb"){
      merged <- mergeWindows(p, tol=5000L, ignore.strand = TRUE)
    }else{
      merged <- mergeWindows(p, tol=2000L, ignore.strand = TRUE)
    }
    cat(length(merged$region))
    dmeta = as.data.frame(p)[, -c(1:7)]
    dmeta[, 1] = 10^(-dmeta$X.log10.pvalue.); dmeta[, 2] = log2(dmeta[,2])
    colnames(dmeta) = c("PValue", "log.fc", "log10.qvalue", "name")
    tabcom <- combineTests(merged$id,  dmeta, pval.col = 1, fc.col = 2 )
    
    pp = merged$region
    elementMetadata(pp) = tabcom
    #p = reduce(p); p10 = reduce(p10);
    #cat(length(p), length(p10), '\n')
    #pp = reduce(pp, drop.empty.ranges=FALSE, min.gapwidth=1000L, with.revmap=TRUE, with.inframe.attrib=FALSE)
    
    peakAnnots = annotatePeak(pp, TxDb=txdb, tssRegion = c(-3000, 3000))
    #xx = data.frame(peakAnnots)
    print(plotAnnoBar(peakAnnots, title = design.matrix$file.name[n]))
    print(plotDistToTSS(peakAnnots, title = design.matrix$file.name[n]))
    
    annotatedPeak = as.data.frame(peakAnnots)
    ### Use the ChIPpeakAnno (we did NOT use here)
    #annotatedPeak <- annotatePeakInBatch(allPeaksList[[1]][[1]], AnnotationData = genes(txdb), output="both")
    #annotatedPeak <- annotatePeakInBatch(allPeaks[[n]], AnnotationData = promoters(genes(txdb), upstream=0, downstream=0))
    ii = grep('geneId', colnames(annotatedPeak))
    mm = match(annotatedPeak[,ii], entrez2symbol$ENTREZID)
    df = data.frame(annotatedPeak, entrez2symbol[mm, ], stringsAsFactors = FALSE)
    kk = which(!is.na(match(df[,1], chroms)));
    df = df[kk, ]
   
    name.annotation = paste0(outDir, "/", sub('\\.xls$', '', design.matrix$file.name[n]), '_combineTests_geneAssignment.txt')
    write.table(df, file = name.annotation, row.name = FALSE, col.name = TRUE, quote=FALSE, sep='\t')
    
  }
  dev.off()
}

##################################################
##################################################
## Section: Differential Analysis without replicates
##################################################
##################################################
## First step:
## Define a set of peaks (the union of peaks) for each factor
## count reads for those peaks and save the tables
Count.Reads.Peak.Union.per.Factor = FALSE
if(Count.Reads.Peak.Union.per.Factor)
{
  ## bam files
  Path4bams = "DATA/BAMs_All_GCc"
  bam.list = list.files(path=Path4bams, pattern = "*.bam$", full.names = TRUE, ignore.case = TRUE)
  
  ## peak files: either use the pre-defined peaks or direct read peak files using design matrix 
  Use.defined.peaks = FALSE;
  if(Use.defined.peaks){
    Path4peaks = "Results/peaks_GCc_merged_with_geneAssignment"
    peak.list = list.files(path=Path4peaks, pattern = "*merged_macs2_broad_fdr_0.1_peakscombineTests_geneAssignment.txt", full.names = TRUE, ignore.case = TRUE)
  }
  
  ## import peak annoation libraries to keep gene assignment for peak unionn 
  library(AnnotationDbi)
  library("ChIPseeker");
  library("rtracklayer");
  library("GenomicFeatures")
  library("GenomicRanges")
  library(Rsubread)
  #library(csaw) # not used csaw function to merge close peaks
  ## import gene annotation
  #txdb <- makeTxDbFromUCSC(genome="mm10", tablename="refGene")
  #saveDb(txdb, file="refGene_mm10.sqlite")
  txdb <- loadDb("refGene_mm10.sqlite")
  #library(AnnotationHub)
  ah = AnnotationHub()
  mm = ah[["AH49583"]]
  entrez2symbol <- unique(select(mm, keys=keys(mm, keytype="ENTREZID"), columns=c("SYMBOL"), keytype="ENTREZID"))
  
  Chr2save = paste0('chr', c(1:19, 'X', 'Y'));
  chroms = Chr2save;
  
  Path2Save = "DATA/read_counts_peaks_All_GCc/"
  if(!dir.exists(Path2Save)){dir.create(Path2Save)}
  
  factors = unique(design.matrix$factor)
  
  for(n in 2:length(factors)) ## loop for each factor
  {
    # n = 1
    cat(factors[n], "\n")
    ## define the union of peaks per factor
    kk = which(design.matrix$factor==factors[n])
    if(length(kk)>1)
    {
      peaks = NULL;
      for(index in 1:length(kk)) 
      {
        k = kk[index];
        p = readPeakFile(design.matrix$file.path[k], as = "GRanges");
        cat(length(p), '\n')
        
        if(factors[n] == "H2AK119"){
          #merged <- mergeWindows(p, tol=5000L, ignore.strand = TRUE)
          merged = reduce(p, drop.empty.ranges=TRUE, min.gapwidth=5000L, with.revmap=TRUE, with.inframe.attrib=FALSE)
        }else{
          #merged <- mergeWindows(p, tol=2000L, ignore.strand = TRUE)
          merged = reduce(p, drop.empty.ranges=TRUE, min.gapwidth=2000L, with.revmap=TRUE, with.inframe.attrib=FALSE)
        }
        cat(length(merged), "\n")
        
        if(index==1){
          peaks = merged; 
        }else{
          peaks= GenomicRanges::union(peaks, merged, ignore.strand=TRUE)
          #peaks10= union(peaks10, p10, ignore.strand=TRUE);
        }
      }
      cat(length(peaks), '\n')
      #mm = match(seqnames(peaks), chroms)
      #xx = peaks[which(!is.na(mm))]
      #peaks = xx;
      #peaks = mergeWindows(peaks, tol = 0L, ignore.strand = TRUE);
      #cat(length(peaks$region), '\n')
      peakAnnots = annotatePeak(peaks, TxDb=txdb, tssRegion = c(-3000, 3000))
      
      ## find genes associated
      annotatedPeak = as.data.frame(peakAnnots)
      ii = grep('geneId', colnames(annotatedPeak))
      mm = match(annotatedPeak[,ii], entrez2symbol$ENTREZID)
      df = data.frame(annotatedPeak, entrez2symbol[mm, ], stringsAsFactors = FALSE)
      df = df[which(!is.na(match(df[,1], chroms))), ]
      
      ## make SAF format annotation for featureCounts
      df$GeneID = paste0(df$seqnames, "_", df$start, "_", df$end)
      SAF = data.frame(GeneID=df$GeneID, Chr=df$seqnames, Start=df$start, End=df$end, Strand="+", assigned.gene=df$SYMBOL, stringsAsFactors = FALSE)
      bams = bam.list[grep(factors[n], bam.list)]
      
      ## count reads for those peak regions
      fc <- featureCounts(files=bams, annot.ext = SAF, countMultiMappingReads = FALSE, minMQS = 10, 
                          ignoreDup = TRUE, strandSpecific = 0, juncCounts = FALSE, nthreads = 8)
      stat = fc$stat;
      counts = fc$counts;
      counts.annot = fc$annotation
      
      res = data.frame(SAF, Length=counts.annot$Length[match(SAF$GeneID, counts.annot$GeneID)], 
                       counts[match(SAF$GeneID, rownames(counts)), ], stringsAsFactors = FALSE)
      
      ## save the read count matrix and also stat for the normalization
      write.table(stat, file = paste0(Path2Save, version.analysis, "_",  factors[n], "_cntStat_chipseq_used4normalization.txt"), 
                  sep='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)
      
      write.table(res, file = paste0(Path2Save, version.analysis, "_",  factors[n], "_readcounts_chipeaks.txt"), 
                  sep='\t', row.names = FALSE, col.names = TRUE, quote=FALSE)
      
      #pp = peaks
      #df <- data.frame(seqnames=seqnames(pp), starts=start(pp)-1, ends=end(pp), names=c(rep(".", length(pp))),
      #                 scores=c(rep(".", length(pp))), strands=strand(pp))
      #kk = which(!is.na(match(df[,1], chroms)));
      #df = df[kk, ]
      #bedname = paste0(path.peaks, nn, '_pval_10_5_PEAK_Union.bed')
      #if(merge.close.peaks){
      #  bedname = paste0(path.peaks, nn, '_pval_10_5_Merge_Close_PEAK_Union.bed')
      #}
      #write.table(df, file=bedname, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    } 
  }
  
}

### Second Step:
### Differential Binding analysis using DEseq2
DB.analysis.using.DEseq2 = FALSE
if(DB.analysis.using.DEseq2)
{
  require('DESeq2')
  path = "DATA/read_counts_peaks_All_GCc"
  xlist <-list.files(path=path, pattern = "*readcounts_chipeaks.txt", full.names = TRUE)
  ylist = list.files(path=path, pattern = "*cntStat_chipseq_used4normalization.txt", full.names = TRUE)
  
  
  ### select factor to analyze
  #factor = 'Ring1B'
  for(factor in factors)
  {
    norms = read.delim(ylist[grep(factor, ylist)], sep='\t', header = TRUE, as.is = 1)
    cts = read.delim(xlist[grep(factor, xlist)], sep='\t', header = TRUE)
    
    #### 
    #### Control the quality 
    raw = as.matrix(cts[, grep("DATA.BAMs", colnames(cts))]);
    rownames(raw) = cts$GeneID;
    raw[which(is.na(raw))] = 0
    mm = match(colnames(raw), colnames(norms))
    norms = norms[, c(1, mm)]
    sample.names = sapply(colnames(raw), function(x) unlist(strsplit(as.character(x), "[.]"))[3], USE.NAMES = FALSE)
    
    colnames(raw) = sample.names;
    colnames(norms)[-1] = sample.names;
    
    #find.samples.conditions(colnames(raw)[1], ID='conditions')
    conds = factor(sapply(colnames(raw), find.samples.conditions, ID='conditions'))
    dds <- DESeqDataSetFromMatrix(raw, DataFrame(conds), ~ conds)
    
    jj = which(norms$Status=="Assigned" | norms$Status=="Unassigned_NoFeatures")
    #ss = as.numeric(norms[match(colnames(dds), norms$sampleName), which(colnames(norms)=="unique.rmdup")]);
    ss = apply(as.matrix(norms[jj, -1]), 2, sum)
    sizeFactors(dds) <- ss/median(ss) 
    #cat(cc, normalizations[nn], "\n")
    #print(sizeFactors(dds))
    
    source("functions_chipSeq.R")
    pdfname = paste0("PLOTs/DB_Analysis_ChIPseq_Qulity_Assessment_", factor , ".pdf")
    pdf(pdfname, width = 12, height = 10)
    
    ## percentages of reads within peaks
    jj = which(norms$Status=="Assigned" | norms$Status=="Unassigned_NoFeatures")
    stat = as.matrix(norms[jj, -1])
    rownames(stat) = c("within.peaks", "others")
    par(mfrow=c(1,1))
    par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(15,3,2,0.2), tcl = -0.3)
    barplot(stat/10^6, horiz = FALSE, names.arg = colnames(stat), las=3, col = c('blue', 'gray'), 
            main='Total nb of reads quantified for features (M)', xlab=NA, legend = rownames(stat), beside=FALSE)
    abline(h=seq(5, 20, by=5), col='red', lty=1, lwd=1.0);#abline(v=c(20, 45), col='red', lty=1, lwd=2.0)
    
    ## Control quality using the RNAseq function
    Assess.ChIPseq.Quality.4DB(dds) 
    
    ## check the replicates
    #toCheck = c("49443", "52740")
    #index = c(); for(cch in toCheck) index = c(index, grep(cch, colnames(fpm)))
    #plot((fpm[, index]+2^-6), log='xy', cex=0.4);
    #abline(0, 1, lwd=2, col='red')
    dev.off()
  
  }
  #########################
  ### pairwise comparison
  dds = estimateDispersions(dds)
  dds = nbinomWaldTest(dds)
  fpm = fpm(dds, robust = TRUE)
  
  
  res1 <- results(dds, contrast = c("conds", "515D10", "AN312"))
  res2 = results(dds, contrast = c("conds", "515D10H3", "AN312"))
  length(which(res1$pvalue<0.05))
  length(which(res2$pvalue<0.05))
  
  plot(res1$baseMean, res1$log2FoldChange, log='x', cex=0.5)
  jj = which(res1$pvalue>0.05); points(res1$baseMean[jj], res1$log2FoldChange[jj], col='blue', cex=0.8)
  
  ### try to remove batch effect
  Batch.consideration = FALSE
  if(Batch.consideration)
  {
    raw = as.matrix(Counts[,-c(1, grep("Input", colnames(Counts)))]);
    rownames(raw) = Counts$peakname;
    raw[which(is.na(raw))] = 0
    require('DESeq2')
    conds = (cbind(rep(c('old','new'), time=3), sapply(colnames(raw), find.samples.conditions, ID='conditions')))
    colnames(conds) = c('batch', 'cell')
    dds <- DESeqDataSetFromMatrix(raw, DataFrame(conds), ~ batch + cell)
    ss = as.numeric(norms[match(colnames(dds), norms$sampleName), which(colnames(norms)=="unique.rmdup")]); 
    sizeFactors(dds) <- ss/median(ss) 
    print(sizeFactors(dds))
    dds = estimateDispersions(dds)
    dds = nbinomWaldTest(dds)
    fpm = fpm(dds, robust = TRUE)
    
    sva.test = FALSE
    if(sva.test){
      library("sva")
      dat  <- counts(dds, normalized = TRUE)
      idx  <- rowMeans(dat) > 1
      dat  <- dat[idx, ]
      mod  <- model.matrix(~conds, colData(dds))
      mod0 <- model.matrix(~   1, colData(dds))
      svseq <- svaseq(dat, mod, mod0=NULL, n.sv = 3)
      
      library(zebrafishRNASeq)
      data(zfGenes)
      filter = apply(zfGenes, 1, function(x) length(x[x>5])>=2)
      filtered = zfGenes[filter,]
      genes = rownames(filtered)[grep("^ENS", rownames(filtered))]
      controls = grepl("^ERCC", rownames(filtered))
      group = as.factor(rep(c("Ctl", "Trt"), each=3))
      dat0 = as.matrix(filtered)
      mod1 = model.matrix(~group)
      mod0 = cbind(mod1[,1])
      svseq = svaseq(dat0,mod1,mod0,n.sv=1)$sv
      plot(svseq,pch=19,col="blue")
    }
  }
  
  
  #rownames(res) = all$gene
  
  if(Filter)
  {
    ccpm = fpm(dds, robust = TRUE)
    mean.2 = apply(ccpm[, c(2, 4)], 1, mean)
    mean.1 = apply(ccpm[, c(1, 3)], 1, mean)
    res = res[which(mean.1>cutoff & mean.2>cutoff), ]
  }
  if(cc=='wt' & nn==1) res0 = res;
  if(cc=='wt' & nn==2) res1 = res;
  
  plot(res$log2FoldChange, -log10(res$pvalue), xlab='log2(FoldChange)', ylab='-log10(pvalue)', cex=0.6, main=paste0(cc, "--", normalizations[nn]))
  abline(v=0, lwd=2.0, col='black')
  abline(h=c(2, 5, 10), lwd=2.0, col='blue')
  text(res$log2FoldChange, -log10(res$pvalue), rownames(res), cex=0.7, offset = 0.3, pos = 3)
  ii = match(c('lsy-6', 'mir-791', 'mir-790', 'mir-793'), rownames(res));
  points(res$log2FoldChange[ii], -log10(res$pvalue)[ii], cex=1.5, col='darkgreen', pch=16)
  
  
  kk0 = grep('common_peaks', bname)
  kk1 = grep('merged', bname);
  kk2 = setdiff(c(1:length(bname)), c(kk0, kk1))
  kk1 = kk1[c(grep(conditions[1], bname[kk1]), grep(conditions[2], bname[kk1]), grep(conditions[3], bname[kk1]))]
  kk2 = kk2[c(grep(conditions[1], bname[kk2]), grep(conditions[2], bname[kk2]), grep(conditions[3], bname[kk2]))]
  xlist = xlist[c(kk0, kk2, kk1)]
  bname = bname[c(kk0, kk2, kk1)]
  
  for(n in 1:nrow(ff)){
    cat(n, '\n')
    xx = readPeakFile(ff$file.path[n], as = "GRanges"); #pps = c(pps, p)
    eval(parse(text = paste0("pp.", n, " = xx")))
  }
  
  TEST.DiffBind = FALSE
  if(TEST.DiffBind)
  {
    BiocStyle::latex()
    savewd <- getwd()
    tmp <-  tempfile(as.character(Sys.getpid()))
    pdf(tmp)
    savewarn <- options("warn")
    options(warn=-1)
    
    library(DiffBind)
    setwd(system.file("extra", package="DiffBind"))
    print(savewd)
    #tamoxifen <- dba(sampleSheet="tamoxifen.csv")
    #tamoxifen <- dba.count(tamoxifen)
    ## tamoxifen <- dba.contrast(tamoxifen)
    ## tamoxifen <- dba.analyze(tamoxifen)
    ## tamoxifen.DB <- dba.report(tamoxifen)
    samples <- read.csv(file.path(system.file("extra", package="DiffBind"),
                                  "tamoxifen.csv"))
    names(samples)
    samples
    ### code chunk number 6: dbaConstruct
    tamoxifen <- dba(sampleSheet="tamoxifen.csv")
    tamoxifen
    plot(tamoxifen)
    tamoxifen <- dba.count(tamoxifen, summits=250)
  }
}
























