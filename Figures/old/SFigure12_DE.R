### libraries
  library(data.table)
  library(ggplot2)
  library(DESeq2)
  library(tidyverse)
  library(gdata)
  library(patchwork)

### load general DE objects
  setwd("/Users/alanbergland/Documents/GitHub/")

  load(file="DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/differential_expression_supplement.Rdata") ### made by `DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/DESeq2/deseq2_QoRTs.R` L129
  load(file="DaphniaPulex20162017Sequencing/AlanFigures/Figure4/differentialExpression_withNames.Rdata") ### made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/Figure4.R`

  #### some data stuff
    ### PCA data
      pcaData <- as.data.table(pcaData)
      pcaData[,clone:=gsub("d8_", "", clone)]
      pcaData[,rep:=tstrsplit(name, "_")[[3]]%>%gsub(".trim.bam", "", .)]

### plot components
  ### PCA
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      theme_bw()

  ### pvalue-distribution
    pval.hist <-
    ggplot(dec, aes(pvalue)) +
    geom_histogram() +
    theme_bw()


layout <- "
AB
CD"

de.mega <-
pcaplot + pval.hist +
cn.plot.lfc + cn.plot.p +
plot_layout(design = layout) +
plot_annotation(tag_levels = 'A')

ggsave(de.mega, file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/SFig12.png")
ggsave(de.mega, file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/SFig12.pdf")






















### copy number DE plot
    cn.plot.lfc <- ggplot() +
    geom_boxplot(data=dec[(cnA-cnB)%in%c(-2, -1, 0, 1, 2)][cnA%in%c(0,1,2) & cnB%in%c(0,1,2)], aes(x=as.factor(I(cnB-cnA)), y=log2FoldChange)) +
    geom_point(data=dec[(cnA-cnB)%in%c(-2, -1, 0, 1, 2)][cnA%in%c(0,1,2) & cnB%in%c(0,1,2)][!is.na(final_QTL_ID)], aes(x=as.factor(I(cnB-cnA)), y=log2FoldChange), color="red") +
    theme_bw() +
    ylab("log2(Fold Change): C / A") +
    xlab("(Copy Number C) - (Copy Number A)")

    cn.plot.p <- ggplot() +
    geom_boxplot(data=dec[(cnA-cnB)%in%c(-2, -1, 0, 1, 2)][cnA%in%c(0,1,2) & cnB%in%c(0,1,2)], aes(x=as.factor(I(cnB-cnA)), y=-log10(pvalue))) +
    geom_point(data=dec[(cnA-cnB)%in%c(-2, -1, 0, 1, 2)][cnA%in%c(0,1,2) & cnB%in%c(0,1,2)][!is.na(final_QTL_ID)], aes(x=as.factor(I(cnB-cnA)), y=-log10(pvalue)), color="red") +
    theme_bw() +
    ylab("-log10(p-value): C / A") +
    xlab("(Copy Number C) - (Copy Number A)")






### mega plot
  pcaplot +
















Daphnia08075
plotCounts(dds, gene=which(rownames(dds)=="Daphnia08075"), returnData=F, intgroup="superclone", transform=T, normalized=T)

### gene expression plots

  gene_counts <- foreach(i=dec[qLFC.noCN.goodChr>.9 & !is.na(old_QTL_ID)]$GeneID, .combine="rbind")%do%{
    gene_counts <- as.data.table(plotCounts(dds, gene=which(rownames(dds)==i), returnData=T, intgroup="superclone", transform=T, normalized=T))
    gene_counts[,gene:=i]
  }
  gene_counts <- merge(gene_counts, dec, by.x="gene", by.y="GeneID")

  table(gene_counts$final_QTL_ID)


  ### QTL 5: 3 genes
  qtl_5_de <-
    ggplot(gene_counts[final_QTL_ID==5],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_wrap(PA42qtl~gene, scales="free") + theme_bw()

  ### QTL 10: 3 genes
  qtl_10_de <-
    ggplot(gene_counts[final_QTL_ID==10],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_grid(PA42qtl~gene, scales="free") + theme_bw()

  ### QTL 12: 4 genes
  qtl_12_de <-
    ggplot(gene_counts[final_QTL_ID==12],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_wrap(PA42qtl~gene, scales="free", nrow=1) + theme_bw()

  ### QTL 7: 1 genes
  qtl_7_de <-
    ggplot(gene_counts[final_QTL_ID==7],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_grid(PA42qtl~gene, scales="free") + theme_bw()

  qtl_9_de <-
    ggplot(gene_counts[final_QTL_ID==9],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_grid(PA42qtl~gene, scales="free") + theme_bw()


qtl_5_de + plot_spacer() + plot_spacer() + plot_spacer() / qtl_8_de / qtl_9_de / qtl_7_de / qtl_9_de


### plot components
  ### PCA
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      theme_bw()

  ###
