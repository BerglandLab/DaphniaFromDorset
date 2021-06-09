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
  library(tidyverse)

### load feature counts
  load("~/featureCounts_QoRTs.Rdata")

### SAF/GTF object
  saf <- flattenGTF("~/Daphnia.aed.0.6.gtf",
      method = "merge")
  saf <- as.data.table(saf)
  setnames(saf, c("Chr", "Start", "End"), c("chr", "start", "end"))
  setkey(saf, GeneID)
  saf <- saf[,list(start=min(start), end=max(end)), list(GeneID, chr)]



### get sample table
  samps <- as.data.frame(fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable"))
  rownames(samps) <- samps$id
  rownames(samps) <- gsub("pe", ".trim", rownames(samps))

### cn.mops output
  load("~/totalresF.Rdata")
  totalresF[,cn.id:=c(1:dim(totalresF)[1])]
  cn <- melt(totalresF, id.vars=c("seqnames", "start", "end", "width", "cn.id"),
              measure.vars=c("April_2017_D8_179_finalmap_mdup.bam",
                             "April_2017_D8_222_finalmap_mdup.bam",
                             "April_2017_D8_349_finalmap_mdup.bam",
                             "May_2017_D8_515_finalmap_mdup.bam"),
              variable.name="clone", value.name="cn")
  cn[,clone:=gsub("_finalmap_mdup.bam", "", clone)]
  cn[,clone:=gsub("April_2017_", "", clone)]
  cn[,clone:=gsub("May_2017_", "", clone)]
  cn[,clone:=tolower(clone)]
  cn <- merge(cn, as.data.table(samps)[,-c("id"),with=F][,list(superclone=superclone[1]), list(clone)], by="clone")
  setnames(cn, c("seqnames"), "chr")

    ### A vs A
    mean(cn[clone=="d8_179"]$cn==cn[clone=="d8_349"]$cn)

    ### C vs C
    mean(cn[clone=="d8_222"]$cn==cn[clone=="d8_515"]$cn)

    ### A vs C
    mean(cn[clone=="d8_179"]$cn==cn[clone=="d8_222"]$cn)
    mean(cn[clone=="d8_349"]$cn==cn[clone=="d8_515"]$cn)

    ggplot(data=cn, aes(x=cn.id, y=clone, fill=cn)) + geom_tile()

  setnames(cn, "cn", "CN")
  cn[,cnum:=gsub("CN", "", CN)%>%as.numeric]

  cn.ag <- cn[,list(cnA=mean(cnum[superclone=="A"]), cnB=mean(cnum[superclone=="C"])),
              list(chr, start, end, cn.id)]

### combine cn.mops output with flattened GTF to flag genes with putative CNVs
  setkey(cn.ag, chr, start, end)
  setkey(saf, chr, start, end)
  cn.gene <- foverlaps(saf, cn.ag, nomatch=NA)
  cn.gene[is.na(cnA), cnA:=2]
  cn.gene[is.na(cnB), cnB:=2]
  cn.gene

  setnames(cn.gene, c("start", "end", "i.start", "i.end"), c("cn_start", "cn_end", "start", "end"))
  setkey(cn.gene, chr, start, end)

### combine cn.gene with QTL peaks file to flag genes under putative QTL
  load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")
  window <- 25000

  peaks.dt <- data.table(qtl=c(1:14), chr=peaks$CHROM, start=peaks$posPeakDeltaSNP-window, end=peaks$posPeakDeltaSNP+window)
  setkey(peaks.dt, chr, start, end)
  cn.gene.qtl <- foverlaps(cn.gene, peaks.dt, nomatch=NA)
  setnames(cn.gene.qtl, c("start", "end", "i.start", "i.end"), c("qtl_start", "qtl_end", "start", "end"))


### deseq
  dds <- DESeqDataSetFromMatrix(countData = fc,
                              colData = samps,
                              design = ~ superclone)


  keep <- rowSums(counts(dds)) >= 50
  table(keep)
  dds <- dds[keep,]

  dds <- DESeq(dds)
  vsd <- vst(dds, blind=F)

  res <- results(dds, format="DataFrame")
  res

  resLFC <- lfcShrink(dds, coef="superclone_C_vs_A", type="apeglm")

  resLFC.dt <- as.data.table(resLFC)
  resLFC.dt[,GeneID:=row.names(resLFC)]

  resLFC.dt[GeneID=="Daphnia00787"]

### combine DE expression w/ CNV data from cn.mops
  de <- merge(resLFC.dt, cn.gene.qtl, by="GeneID", all.x=T)
  de[GeneID=="Daphnia00787"]

### flag genes that are not on the 12 major chromosomes
  de.ag <- de[,.N, chr]
  de.ag[,goodChr:=I(N>100)]
  de <- merge(de, de.ag, by="chr")

  de[,qLFC:=rank(abs(log2FoldChange))/(length(log2FoldChange)+1)]
  de[cnA==2 & cnB==2,qLFC.noCN:=rank(abs(log2FoldChange))/(length(log2FoldChange)+1)]
  de[cnA==2 & cnB==2 & goodChr==T,qLFC.noCN.goodChr:=rank(abs(log2FoldChange))/(length(log2FoldChange)+1)]
  de[cnA==2 & cnB==2 & goodChr==T,qp.noCN.goodChr:=rank(abs(pvalue))/(length(pvalue)+1)]

  de[GeneID=="Daphnia00787"]

### add in PC loading values
  pcaData <- plotPCA(vsd, intgroup=c("superclone", "clone"), returnData=T, ntop=10000)
  t.test(PC1~superclone, pcaData)
  save(pcaData, res, dds, file="DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/differential_expression_supplement.Rdata")


  ggplot(data=pcaData, aes(x=PC1, y=PC2, color=superclone, shape=clone)) + geom_point()

  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(100000,
      length(rv)))]
  v <- assay(vsd)

  res.pca <- prcomp(t(v[row.names(v)%in%de[cnA==2 & cnB==2 & goodChr==T]$GeneID, ]))
  hm <- get_pca_var(res.pca)
  contrib <- hm$contrib
  con.dt <- as.data.table(contrib)
  con.dt[,GeneID:=dimnames(hm$contrib)[[1]]]
  dec <- merge(de, con.dt, all.x=T, by="GeneID")
  dec[GeneID=="Daphnia00787"]

### re-calculate q-values
  dec[cnA==2 & cnB==2 & goodChr==T, pa:=p.adjust(pvalue)]
  dec[cnA==2 & cnB==2 & goodChr==T, qval:=p.adjust(pvalue, method="fdr")]


### top DE-QTL genes
  dec[cnA==2 & cnB==2 & goodChr==T]
  dec[!is.na(qtl)][pa<.05]
  chisq.test(!is.na(dec$qtl), dec$pa<.05)
  table(dec[!is.na(qtl)][qval<.05]$qtl)

### save
  save(dec, file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/differentialExpression.Rdata")






ggplot(data=dec[cnA==2 & cnB==2 & goodChr==T][order(qtl, na.last=F)],
            aes(y=-log10(pvalue), x=log2FoldChange, color=as.factor(!is.na(qtl)))) +
geom_point(alpha=.49)






  ggplot(data=dec[order(qtl, na.last=F)],
        aes(x=qLFC.noCN.goodChr, y=log2FoldChange,
            color=as.factor(qtl))) +
  geom_point() +
  theme_bw() +
  geom_vline(xintercept=.9)

  fisher.test(table(dec$padj<=.05, dec$cnA!=dec$cnB))
