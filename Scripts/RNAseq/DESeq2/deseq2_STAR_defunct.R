### scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/rnaseq/featureCounts_trim.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf ~/.
### scp aob2x@rivanna.hpc.virginia.edu:~/featureCounts_STAR.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:~/genes_rd_missing.Rdata ~/.

### libraries
  library(data.table)
  library(factoextra)

  library(DESeq2)
  library(ggplot2)
  library(Rsubread)
  library("RColorBrewer")
  library(pheatmap)
  library(patchwork)

### load feature counts
  load("~/featureCounts_STAR.Rdata")

### load gene read depth info
  load("~/genes_rd_missing.Rdata")

### SAF/GTF object
  saf <- flattenGTF("~/Daphnia.aed.0.6.gtf",
      method = "merge")

### get sample table
  samps <- as.data.frame(fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable"))
  rownames(samps) <- samps$id
  rownames(samps) <- gsub("pe", ".trim", rownames(samps))
  dds <- DESeqDataSetFromMatrix(countData = fc,
                              colData = samps,
                              design = ~ superclone)


  keep <- rowSums(counts(dds)) >= 10
  table(keep)
  dds <- dds[keep,]

  dds <- DESeq(dds)
  vsd <- vst(dds, blind=FALSE)

  res <- results(dds, format="DataFrame")
  res

  resLFC <- lfcShrink(dds, coef="superclone_C_vs_A", type="apeglm")

### "best QTL gene"
  plotCounts(dds, gene=which(rownames(dds)=="Daphnia11267"), intgroup="superclone")

  gene_counts <- plotCounts(dds, gene=which(rownames(dds)=="Daphnia00787"), returnData=T, intgroup="superclone", transform=T, normalized=T)
  geneCounts <- ggplot(gene_counts, aes(x=superclone, y=count)) + geom_point() + ylab("Normalized expression - logscale")

### volcano, and LFC density plot
  resLFC$GeneId <- rownames(resLFC)
  res.dt <- as.data.table(resLFC)

  res.dt[GeneId=="Daphnia00787"]
  mean(res.dt$pvalue <= res.dt[GeneId=="Daphnia00787"]$pvalue, na.rm=T)
  mean(abs(res.dt$log2FoldChange) >= abs(res.dt[GeneId=="Daphnia00787"]$log2FoldChange), na.rm=T)

### PCA analysis
  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(100000,
      length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]))
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
        name = colnames(vsd))
  pcaData <- plotPCA(vsd, intgroup=c("superclone", "clone"), returnData=F, ntop=10000)

  plot(pca$rotation[which(dimnames(pca$rotation)[1]=="Daphnia00787"),1] ~ pca$rotation[,2])


### combine DE expression w/ read-depth and missing data rates
  genes.ag.ag <- genes.ag[,list(scA=mean(dp.norm.mu[superclone=="A"], na.rm=T) , scB=mean(dp.norm.mu[superclone=="C"], na.rm=T),
                                  missingA=mean(missing.rate[superclone=="A"], na.rm=T) , missingB=mean(missing.rate[superclone=="C"], na.rm=T)), list(GeneId=GeneID)]
  res.gene <- merge(res.dt, genes.ag.ag, by="GeneId")
  res.gene[,delta:=scA-scB]
  res.gene[,quanLFC:=rank(abs(log2FoldChange))/(length(log2FoldChange)-1)]

### combine with QTL peaks
### combine with QTL
  load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")
  peaks.dt <- data.table(qtl=c(1:14), Chr=peaks$CHROM, Start=peaks$posPeakDeltaSNP-15000, End=peaks$posPeakDeltaSNP+15000)

  saf.dt <- as.data.table(saf)
  setkey(saf.dt, GeneID)
  saf.dt<-unique(saf.dt, by = key(saf.dt), fromLast = TRUE)

  setkey(peaks.dt, Chr, Start, End)
  setkey(saf.dt, Chr, Start, End)

  qtl.genes <- foverlaps(peaks.dt, saf.dt)

  setnames(qtl.genes, "GeneID", "GeneId")
  setkey(qtl.genes, "GeneId")

  resLFC$GeneId <- row.names(resLFC)
  resLFC <- as.data.table(resLFC)

  res.qtl <- merge(resLFC, qtl.genes[,c("GeneId", "qtl"),with=F], by="GeneId", all.x=T)
  setkey(res.qtl, "GeneId")

  res.qtl[!is.na(qtl)][order(padj)]
  setnames(saf.dt, "GeneID", "GeneId")
  res.qtl <- merge(res.qtl, saf.dt, by="GeneId")

  res.qtl.ag <- res.qtl[,list(.N), Chr]

  dim(res.qtl)
  res.qtl <- res.qtl[Chr%in%res.qtl.ag[1:12]$Chr]
  dim(res.qtl)

  res.qtl[,mid:=Start/2+End/2]
  res.qtl[,quan:=rank(padj)/(length(padj)-1)]

  res.qtl.gene <- merge(res.qtl, genes.ag.ag, by="GeneId", all.x=T)
  res.qtl.gene <- res.qtl.gene[missingA<1 & missingB<1]
  res.qtl.gene[,quan:=rank(padj)/(length(padj)-1)]



### save
  save(res.gene, vsd, pcaData, res.dt, gene_counts, dds, genes.ag, res.qtl, file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/STAR/DESeq_star.output")






### #OLD
  ### volcano combined with genes.ag
    genes.ag.ag <- genes.ag[,list(scA=mean(dp.norm.mu[superclone=="A"], na.rm=T) , scB=mean(dp.norm.mu[superclone=="C"], na.rm=T),
                                    missingA=mean(missing.rate[superclone=="A"], na.rm=T) , missingB=mean(missing.rate[superclone=="C"], na.rm=T)), list(GeneId=GeneID)]
    res.gene <- merge(res.dt, genes.ag.ag, by="GeneId")
    res.gene[,delta:=scA-scB]
    res.gene[,quanLFC:=rank(abs(log2FoldChange))/(length(log2FoldChange)-1)]

    fisher.test(table(res.gene$quanLFC>.995, res.gene$missingA==1 | res.gene$missingB==1))

    ggplot(data=res.gene, aes(y=delta, x=log2FoldChange)) + geom_point()

    ggplot(data=res.gene, aes(x=scA, y=scB, color=log2FoldChange)) + geom_point()


    fisher.test(table(res.gene$missingB==1 & res.gene$missingA<1, res.gene$log2FoldChange< -4))

    summary(lm(log2FoldChange~scA*scB, res.gene))



### combined plot
  layout <- "
  AABB
  CDEE
  "

  pcaplot + volcano + magicgene_lfc + magicgene_basemean + geneCounts +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A') +
  ggtitle("Daphnia00787")


### combine with QTL
  load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")
  peaks.dt <- data.table(qtl=c(1:14), Chr=peaks$CHROM, Start=peaks$posPeakDeltaSNP-50000, End=peaks$posPeakDeltaSNP+50000)

  saf.dt <- as.data.table(saf)
  setkey(saf.dt, GeneID)
  saf.dt<-unique(saf.dt, by = key(saf.dt), fromLast = TRUE)

  setkey(peaks.dt, Chr, Start, End)
  setkey(saf.dt, Chr, Start, End)

  qtl.genes <- foverlaps(peaks.dt, saf.dt)

  setnames(qtl.genes, "GeneID", "GeneId")
  setkey(qtl.genes, "GeneId")

  resLFC$GeneId <- row.names(resLFC)
  resLFC <- as.data.table(resLFC)

  res.qtl <- merge(resLFC, qtl.genes[,c("GeneId", "qtl"),with=F], by="GeneId", all.x=T)
  setkey(res.qtl, "GeneId")

  res.qtl[!is.na(qtl)][order(padj)]
  setnames(saf.dt, "GeneID", "GeneId")
  res.qtl <- merge(res.qtl, saf.dt, by="GeneId")

  res.qtl.ag <- res.qtl[,list(.N), Chr]

  dim(res.qtl)
  res.qtl <- res.qtl[Chr%in%res.qtl.ag[1:12]$Chr]
  dim(res.qtl)

  res.qtl[,mid:=Start/2+End/2]
  res.qtl[,quan:=rank(padj)/(length(padj)-1)]

  res.qtl.gene <- merge(res.qtl, genes.ag.ag, by="GeneId", all.x=T)
  res.qtl.gene <- res.qtl.gene[missingA<1 & missingB<1]
  res.qtl.gene[,quan:=rank(padj)/(length(padj)-1)]



  res.qtl.ag <- res.qtl.gene[,list(frac=mean(quan<=.1),
                              TT=sum(quan<=.1),
                              expected=.1*.N,
                              .N,
                              p=binom.test(sum(quan<=.1), sum(quan>.1), .1)$p.value), list(qtl)]

  res.qtl.gene[qtl==10][quan<=.1]

### plot
  res.qtl[,mid:=Start/2+End/2]
  res.qtl[,quan:=rank(padj)/(length(padj)-1)]

  ggplot(data=res.qtl[order(quan, na.last=F)], aes(x=abs(log2FoldChange), y=pvalue, color=!is.na(qtl),
        shape=padj<1e-20)) +
  geom_point(size=1, alpha=.75)


  ggplot(data=res.qtl[order(qtl, na.last=F)], aes(x=log2FoldChange, y=-log10(pvalue), color=!is.na(qtl))) + geom_point()


    ggplot(data=res.qtl[order(qtl, na.last=F)], aes(x=mid, y=-log10(padj), color=!is.na(qtl),
          shape=padj<1e-10)) +
    geom_point(size=.75, alpha=.75) + facet_wrap(~Chr, nrow=1, scales="free_x")



  scp aob2x@rivanna.hpc.virginia.edu:~/res_deseq.Rdata ~/.


  library(data.table)
  library(ggplot2)
  load("~/res_deseq.Rdata")

  ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj))) + geom_point()
