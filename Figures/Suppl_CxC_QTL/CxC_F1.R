### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(patchwork)
  library(lme4)


### setwd
  setwd("/Users/alanbergland/Documents/GitHub")

### load F1 cross data (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")

### load in PoolSeq data (made by `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLSeqR.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")

### load in F1 mapping data
  #load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")
  #load("~/lme4qtl_output.AxC.long.Rdata")
  load("DaphniaPulex20162017Sequencing/AlanFigures/Suppl_CxC_QTL/lme4qtl_output.CxC.long.SUMMARIZED.Rdata")


### Load in PA42 chr key
  PA42 <- fread("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/5kbchrassignHiCnew.csv")
  PA42$hybridchr <- paste("Chr", PA42$PA42chr, "\n", paste("Scaff", "\n", tstrsplit(PA42$chr, "_")[[2]], sep=""), sep="")
  PA42sub <- PA42[, c("chr", "PA42chr", "hybridchr")]

### some PoolSeq stuff
  setnames(peaks, "CHROM", "chr")
  setnames(gprime, "CHROM", "chr")
  peaks[,old_QTL_ID:=c(1:14)]

  setkey(PA42sub, chr)
  setkey(peaks, chr)
  setkey(gprime, chr)
  mpeaks <- merge(peaks, PA42sub)
  mgprime <- merge(gprime, PA42sub)

  mgprime[,Gprime.y:=Gprime*ifelse(rep==1, 1, -1)]
  mpeaks[,maxGprime.y:=maxGprime*ifelse(rep==1, 1, -1)]

  mgprime[,pos:=POS]
  setkey(mgprime, hybridchr, pos)
  merge(mgprime,
        data.table(hybridchr=mpeaks$hybridchr,
                   pos=mpeaks$posMaxGprime, key="hybridchr,pos"))[,]

  mgprime$hybridchr <- factor(mgprime$hybridchr, levels=c("Chr1\nScaff\n1931", "Chr2\nScaff\n9198",
    "Chr3\nScaff\n9199", "Chr4\nScaff\n9197", "Chr5\nScaff\n9200", "Chr6\nScaff\n2373",
    "Chr7\nScaff\n7757", "Chr8\nScaff\n6786", "Chr9\nScaff\n1863", "Chr10\nScaff\n2217",
    "Chr11\nScaff\n9201", "Chr12\nScaff\n2158"))

  mpeaks$hybridchr <- factor(mpeaks$hybridchr, levels=c("Chr1\nScaff\n1931", "Chr2\nScaff\n9198",
    "Chr3\nScaff\n9199", "Chr4\nScaff\n9197", "Chr5\nScaff\n9200", "Chr6\nScaff\n2373",
    "Chr7\nScaff\n7757", "Chr8\nScaff\n6786", "Chr9\nScaff\n1863", "Chr10\nScaff\n2217",
    "Chr11\nScaff\n9201", "Chr12\nScaff\n2158"))


  setkey(mgprime, hybridchr)
  setkey(mpeaks, rep, hybridchr)

  setkey(mpeaks, PA42chr, posMaxGprime)
  mpeaks[,final_QTL_ID:=c(1:14)]


##########################
### F1 mapping results ###
##########################
  setkey(o.ag.plot, chr)
  setkey(PA42, chr)
  mo.ag.plot <- merge(o.ag.plot, PA42)

  o.ag.perm <- mo.ag.plot[,list(p.aov.thr=mean(p.aov.thr)), list(term, chr, hybridchr, q)]

  o.ag.perm$hybridchr <- factor(o.ag.perm$hybridchr, levels=c("Chr1\nScaff\n1931", "Chr2\nScaff\n9198",
    "Chr3\nScaff\n9199", "Chr4\nScaff\n9197", "Chr5\nScaff\n9200", "Chr6\nScaff\n2373",
    "Chr7\nScaff\n7757", "Chr8\nScaff\n6786", "Chr9\nScaff\n1863", "Chr10\nScaff\n2217",
    "Chr11\nScaff\n9201", "Chr12\nScaff\n2158"))

  mo.ag.plot$hybridchr <- factor(mo.ag.plot$hybridchr, levels=c("Chr1\nScaff\n1931", "Chr2\nScaff\n9198",
    "Chr3\nScaff\n9199", "Chr4\nScaff\n9197", "Chr5\nScaff\n9200", "Chr6\nScaff\n2373",
    "Chr7\nScaff\n7757", "Chr8\nScaff\n6786", "Chr9\nScaff\n1863", "Chr10\nScaff\n2217",
    "Chr11\nScaff\n9201", "Chr12\nScaff\n2158"))


#### plo
    f1.plot <-
    ggplot() +
    geom_vline(data=mpeaks, aes(xintercept=posMaxGprime, linetype=as.factor(rep))) +
    geom_line(data=mo.ag.plot, aes(x=pos, y=-log10(p.aov), color=hybridchr), size=1.5) +
    geom_point(data=mo.ag.plot[p.aov<=p.aov.thr], aes(x=pos, y=-log10(p.aov)), size=.5) +
    geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.aov.thr), linetype=as.factor(q))) +
    facet_grid(term~hybridchr, scales="free_x", space = "free_x") +
    theme_bw() + theme(strip.text.x = element_text(size = 8)) +
    scale_x_continuous(breaks=seq(from=1, to=13, by=2)*10^6, labels=seq(from=1, to=13, by=2)) +
    theme(legend.position = "none") + xlab("Position (Mb)") + ylab("-log10(p)")


  ### save
    ggsave(f1.plot, file="DaphniaPulex20162017Sequencing/AlanFigures/Suppl_CxC_QTL/CxC_f1.pdf")
    ggsave(f1.plot, file="DaphniaPulex20162017Sequencing/AlanFigures/Suppl_CxC_QTL/CxC_f1.png")
