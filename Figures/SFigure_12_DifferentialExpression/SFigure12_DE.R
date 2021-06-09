### libraries
  library(data.table)
  library(ggplot2)
  library(DESeq2)
  library(tidyverse)
  library(gdata)
  library(patchwork)
  library(foreach)

### load general DE objects
  setwd("/Users/alanbergland/Documents/GitHub/")

  load(file="DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/differential_expression_supplement.Rdata") ### made by `DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/DESeq2/deseq2_QoRTs.R` L129
  load(file="DaphniaPulex20162017Sequencing/AlanFigures/Figure4/differentialExpression_withNames.Rdata") ### made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/Figure4.R`

  #### some data stuff
    ### PCA data
      pcaData <- as.data.table(pcaData)
      pcaData[,clone:=gsub("d8_", "", clone)]
      pcaData[,rep:=tstrsplit(name, "_")[[3]]%>%gsub(".trim.bam", "", .)]

    ### qtl plots
      gene_counts <- foreach(i=dec[qLFC.noCN.goodChr>.9 & !is.na(old_QTL_ID)]$GeneID, .combine="rbind")%do%{
        gene_counts <- as.data.table(plotCounts(dds, gene=which(rownames(dds)==i), returnData=T, intgroup="superclone", transform=T, normalized=T))
        gene_counts[,gene:=i]
      }
      gene_counts <- merge(gene_counts, dec, by.x="gene", by.y="GeneID")

      table(gene_counts$final_QTL_ID)


### plot components
  ### PCA
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    pcaplot <- ggplot(pcaData, aes(PC1, PC2, color=superclone, shape=clone)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      theme_bw() +
      theme(strip.text = element_text(size = 14),
            axis.text=element_text(size=12),
            axis.title=element_text(size=14))


  ### pvalue-distribution
    pval.hist <-
    ggplot(dec, aes(pvalue)) +
    geom_histogram() +
    theme_bw() +
    theme(strip.text = element_text(size = 14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))

  ### gene expression plots
    qtl_12_de <-
    ggplot(gene_counts[final_QTL_ID==12],
    aes(x=superclone, y=count)) +
    geom_point() +
    ylab("Normalized expression - logscale") +
    facet_wrap(~gene, scales="free", nrow=1) +
    theme_bw() +
    theme(strip.text = element_text(size = 14),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))


### mega plot
  layout <- "
  AB
  CC"

  de.mega <-
  pcaplot + pval.hist +
  qtl_12_de +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')

  de.mega

  ggsave(de.mega, file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/SFig12.png")
  ggsave(de.mega, file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/SFig12.pdf")






















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
