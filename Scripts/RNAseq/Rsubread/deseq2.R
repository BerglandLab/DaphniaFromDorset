### scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/rnaseq/featureCounts_trim.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf ~/.

### libraries
  library(data.table)
  library(DESeq2)
  library(ggplot2)
  library(Rsubread)
  library("RColorBrewer")
  library(pheatmap)
  library(patchwork)

### load feature counts
  load("~/featureCounts_trim.Rdata")

### SAF/GTF object
  saf <- flattenGTF("~/Daphnia.aed.0.6.gtf",
      method = "merge")

### get sample table
  samps <- as.data.frame(fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable"))
  rownames(samps) <- samps$id
  rownames(samps) <- gsub("pe", ".trim", rownames(samps))
  dds <- DESeqDataSetFromMatrix(countData = fc$counts,
                              colData = samps,
                              design = ~ superclone)


  keep <- rowSums(counts(dds)) >= 50
  table(keep)
  dds <- dds[keep,]

  dds <- DESeq(dds)
  #vsd <- vst(dds, blind=FALSE)

  res <- results(dds, format="DataFrame")
  res

  resLFC <- lfcShrink(dds, coef="superclone_C_vs_A", type="apeglm")

### "best QTL gene"
  plotCounts(dds, gene=which(rownames(dds)=="Daphnia00796"), intgroup="superclone")
  gene_counts <- plotCounts(dds, gene=which(rownames(dds)=="Daphnia00787"), returnData=T, intgroup="superclone", transform=T, normalized=T)
  geneCounts <- ggplot(gene_counts, aes(x=superclone, y=count)) + geom_point() + ylab("Normalized expression - logscale")


  resLFC$GeneId <- rownames(resLFC)
  res.dt <- as.data.table(resLFC)

  res.dt[GeneId=="Daphnia00787"]
  mean(res.dt$pvalue <= res.dt[GeneId=="Daphnia00787"]$pvalue, na.rm=T)
  mean(abs(res.dt$log2FoldChange) >= abs(res.dt[GeneId=="Daphnia00787"]$log2FoldChange), na.rm=T)

  volcano <- ggplot() +
  geom_point(data=res.dt, aes(x=log2FoldChange, y=-log10(pvalue)), color="red") +
  geom_point(data=res.dt[GeneId=="Daphnia00787"], aes(x=log2FoldChange, y=-log10(pvalue)), color="blue")

  magicgene <- ggplot(data=res.dt, aes(abs(log2FoldChange))) +
  geom_density() +
  geom_vline(xintercept=abs(res.dt[GeneId=="Daphnia00787"]$log2FoldChange), color="red") +
  geom_text(data=data.frame(x=4, y=1, label=mean(abs(res.dt$log2FoldChange) >= abs(res.dt[GeneId=="Daphnia00787"]$log2FoldChange), na.rm=T)),
            aes(x=x, y=y,
                label=round(label, 3)))


### PCA analysis
  pcaData <- plotPCA(vsd, intgroup=c("superclone", "clone"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))



### clustering
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  heatmap <- pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)


### combined plot
  layout <- "
  AABB
  CCDD
  "

  pcaplot + volcano + magicgene + geneCounts +
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

### plot
  res.qtl[,mid:=Start/2+End/2]

  ggplot(data=res.qtl[padj<1e-10][order(qtl, na.last=F)], aes(x=mid, y=log2FoldChange, color=!is.na(qtl),
        shape=padj<1e-20)) +
  geom_point(size=1, alpha=.75) + facet_wrap(~Chr, nrow=1, scales="free_x")


  ggplot(data=res.qtl[order(qtl, na.last=F)], aes(x=log2FoldChange, y=-log10(pvalue), color=!is.na(qtl))) + geom_point()


    ggplot(data=res.qtl[order(qtl, na.last=F)], aes(x=mid, y=-log10(padj), color=!is.na(qtl),
          shape=padj<1e-10)) +
    geom_point(size=.75, alpha=.75) + facet_wrap(~Chr, nrow=1, scales="free_x")



  scp aob2x@rivanna.hpc.virginia.edu:~/res_deseq.Rdata ~/.


  library(data.table)
  library(ggplot2)
  load("~/res_deseq.Rdata")

  ggplot(data=res, aes(x=log2FoldChange, y=-log10(padj))) + geom_point()
