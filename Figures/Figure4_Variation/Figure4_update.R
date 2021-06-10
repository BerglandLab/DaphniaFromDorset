### libraries
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(ggrepel)

############
### Data ###
############

### setwd
  #setwd("/Users/kbkubow/Documents/GitHub")
  setwd("/Users/alanbergland/Documents/GitHub")

### load F1 cross data (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")

### load in PoolSeq data (made by `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLSeqR.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")

### load in F1 mapping data
  #load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")
  #load("~/lme4qtl_output.AxC.long.Rdata")
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.SUMMARIZED.Rdata")

### Load in PA42 chr key
  PA42 <- fread("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/5kbchrassignHiCnew.csv")

### Load in qtl_polar (wild daps only: `DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_superclone_test/analysis.polarizeHapolotypes_poolSeq.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/qtl_polar.Rdata")

### load in F1 validation (`DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_superclone_test/F1_mapping_validation.polarizeHapolotypes_poolSeq.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/f1_pool_polar.Rdata")
  f1.val <- f1.pool.merge[,list(propmale=mean(propmale), N=mean(N),
                                nMale=sum(geno=="male_pe") + 2*sum(geno=="male_male")),
                           list(clone, gr)]

  f1.val.glm <- glm(propmale~nMale, f1.val, weights=N, family="binomial")
  with(summary(f1.val.glm), 1 - deviance/null.deviance)

### load basic expression data (`DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/DESeq2/deseq2_QoRTs.R`)
  load(file="DaphniaPulex20162017Sequencing/AlanFigures/Figure4/differentialExpression.Rdata")

### Load GEVA output (`DaphniaPulex20162017Sequencing/AlanFigures/Figure4/GEVA_figure/GEVA_plot.R`)
  # scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/age.tar.gz ~/.
  load(file="DaphniaPulex20162017Sequencing/AlanFigures/Figure4/tmrca.output")

############################
### male proportion plot ###
############################

  male$gr <- ifelse(male$SC=="selfedA", "A", ifelse(male$SC=="B", "B", male$gr))
  male.ag <- male[!is.na(gr),list(propmale=sum(Males)/sum(NewTotal),
                                  N=sum(NewTotal)),
                    list(clone, gr)]
  male.ag[,se:=sqrt((propmale*(1-propmale))/N)]

  male.ag[,lci:=propmale-1.96*se]
  male.ag[,uci:=propmale+1.96*se]

  male.ag$gr <- factor(male.ag$gr, levels=c("A", "AxC", "C", "CxC"))


### Fill Rate
  epp$gr <- ifelse(epp$SC=="selfedA", "A", ifelse(epp$SC=="B", "B", epp$gr))
  epp$counts <- c(1)

  epp.ag <- epp[!is.na(gr),list(fillrate=sum(fill*TotalEppB, na.rm=T)/sum(TotalEppB),
                                N=sum(TotalEppB), meanEMB=mean(TotalEmb),
                                sdEMB=sd(TotalEmb), sampsize=sum(counts),
                                meanEpp=mean(TotalEpp), sdEpp=sd(TotalEpp)),
                  list(clone, gr)]

  epp.ag[,se:=sqrt((fillrate*(1-fillrate))/N)]
  epp.ag[,lci:=fillrate-1.96*se]
  epp.ag[,uci:=fillrate+1.96*se]

  epp.ag$gr <- factor(epp.ag$gr, levels=c("A", "AxC", "C", "CxC"))


#########################
### Pool seq results. ###
#########################
    setnames(peaks, "CHROM", "chr")
    setnames(gprime, "CHROM", "chr")
    peaks[,old_QTL_ID:=c(1:14)]

    PA42$hybridchr <- paste("Chr", PA42$PA42chr, "\n", paste("Scaff", "\n", tstrsplit(PA42$chr, "_")[[2]], sep=""), sep="")
    PA42sub <- PA42[, c("chr", "PA42chr", "hybridchr")]
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


#############################
### QTL polarization plot ###
#############################
  qtlcoding <- data.table(qtl=c(1:14), PA42qtl=c(3, 4, 2, 1, 9, 10, 11, 13, 14, 12, 8, 6, 7, 5))

  mqtl.polar <- merge(qtl.polar, qtlcoding, by="qtl")
  qtl.polar <- mqtl.polar
  setkey(mqtl.polar, PA42qtl)

  qtl.polar$geno <- ifelse(qtl.polar$n.male_male > qtl.polar$n.male_pe &
    qtl.polar$n.male_male > qtl.polar$n.pe_pe, "male+/male+",
    ifelse(qtl.polar$n.male_pe > qtl.polar$n.male_male &
    qtl.polar$n.male_pe > qtl.polar$n.pe_pe, "male+/male-",
    ifelse(qtl.polar$n.pe_pe > qtl.polar$n.male_male &
    qtl.polar$n.pe_pe > qtl.polar$n.male_pe,"male-/male-", "other")))

  qtlsub <- qtl.polar[SC.uniq=="C" & year=="2017" & pond.y=="D8" | SC.uniq=="A" &
    year=="2017" & pond.y=="D8"| SC.uniq=="B" & year=="2018"]

  mq <- merge(male.ag, qtl.polar, by="clone")


#######################
### gene expression ###
#######################

  dec <- merge(dec, qtlcoding, by="qtl", all.x=T)
  setnames(dec, "qtl", "old_QTL_ID")
  dec <- merge(dec, mpeaks[,c("chr", "posPeakDeltaSNP", "PA42chr", "hybridchr", "final_QTL_ID", "old_QTL_ID"), with=F], all.x=T, by="old_QTL_ID")

  save(dec, file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/differentialExpression_withNames.Rdata")

  dec[,lab:=""]
  dec[qLFC.noCN.goodChr>.9 & !is.na(old_QTL_ID), lab:=paste("QTL-", final_QTL_ID, ": ", GeneID, sep="")]
  dec[qLFC.noCN.goodChr>.9 & !is.na(old_QTL_ID), lab2:=final_QTL_ID]

  table(dec[lab!=""]$PA42qtl)


########################
### tmrca expression ###
########################
  tmrca[,chr:=CHROM.x]
  setkey(mpeaks, chr, posPeakDeltaSNP)
  setkey(tmrca, chr, posPeakDeltaSNP)

  tmrca <- merge(tmrca, mpeaks[,c("chr", "posPeakDeltaSNP", "PA42chr", "hybridchr", "final_QTL_ID"), with=F], all.x=T)
  tmrca[,lab:=""]
  tmrca[!is.na(hybridchr), lab:=final_QTL_ID]


#########################
### figure components ###
#########################

  male.plot <- ggplot(data=male.ag[!is.na(gr)], aes(x=gr, y=propmale)) +
  geom_linerange(aes(ymin=lci, ymax=uci), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Male proportion") +theme_bw() + theme(legend.position="none")


  epp.plot <- ggplot(data=epp.ag[!is.na(gr)], aes(x=gr, y=fillrate)) +
  geom_linerange(aes(ymin=lci, ymax=uci), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Ephippial fill rate") +
  theme_bw() + theme(legend.position="none")

  poolseq <- ggplot() +
  #geom_vline(data=mpeaks, aes(xintercept=posMaxGprime), color="black") +
  geom_line(data=mgprime, aes(x=pos, y=Gprime.y, group=rep, color=as.factor(rep)), size=1.5) +
  geom_segment(data=mpeaks, aes(x=posMaxGprime, xend=posMaxGprime,
                            y=maxGprime.y+0.15*ifelse(rep==1, 1, -1), yend=5*ifelse(rep==1, 1, -1))) +
  geom_hline(data=mgprime[,list(minG=min(Gprime[qvalue<=0.05])), list(rep)],
         aes(yintercept=minG*ifelse(rep==1, 1, -1))) +
  geom_text(data=mpeaks, aes(x=posMaxGprime, y=6*ifelse(rep==1, 1, -1), label=final_QTL_ID), size=3) +
  facet_grid(~hybridchr, scales="free_x", space = "free_x") +
  theme(legend.position = "none") +
  theme_bw() +
  scale_x_continuous(breaks=seq(from=1, to=13, by=2)*10^6, labels=seq(from=1, to=13, by=2)) +
  scale_y_continuous(breaks=c(-4, 0, 4), labels=c(4,0,4), limits=c(-6.5, 6.5)) +
  theme(strip.text.x = element_text(size = 8), legend.position="none") +
  ylab("Gprime") + xlab("Position (Mb)") +
  labs(color="Replicate") +
  ylab("G'")


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


  qtlpolar.strip <-
  ggplot(data=qtlsub[SC.uniq!="B"], aes(x=factor(SC.uniq, levels=c("C", "A")), y=PA42qtl, fill=geno)) +
  geom_tile(color="black",size=.5, width=0.95, height=0.95) +
  coord_flip() +
  xlab("Clonal lineage") +
  ylab("QTL") +
  labs(fill="Genotype") +
  theme_bw() +
  scale_y_continuous(breaks=c(1:14), expand=c(0,.15)) +
  theme(legend.position="bottom")


  f1polar.plot <-
  ggplot() +
  geom_point(data=f1.val, aes(x=nMale, y=propmale, color=factor(gr, levels=c("A", "AxC", "C", "CxC")))) +
  theme_bw() +
  xlab("Num. male alleles") +
  ylab("Proportion.male") +
  geom_text(aes(x=14, y=.15, label=paste("z =", round(summary(f1.val.glm)$coef[2,3], 1))), size=2) +
  geom_text(aes(x=14, y=.13, label=paste("p =", 5*10^-26)), size=3.5) +
  theme(legend.position="none")


  geneexpression.plot <-
  ggplot(data=dec[cnA==2 & cnB==2 & goodChr==T][order(final_QTL_ID, na.last=F)],
            aes(y=-log10(pvalue), x=log2FoldChange, color=as.factor(!is.na(final_QTL_ID)), label=lab2)) +
  geom_point(alpha=.49) +
  theme(legend.position = "none") +
  theme_bw() +
  geom_label_repel(
    data=dec[cnA==2 & cnB==2 & goodChr==T][order(final_QTL_ID, na.last=F)][PA42qtl==8],
    force   = 10, # do not pull toward data points
    nudge_y      = 75,
    nudge_x       =5,
    direction    = "y",
    angle        = 0,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    size=3.5, color="black"
  ) +
  geom_label_repel(
    data=dec[cnA==2 & cnB==2 & goodChr==T][order(final_QTL_ID, na.last=F)][PA42qtl==10],
    force   = 10, # do not pull toward data points
    nudge_y      = 100,
    nudge_x       =0,
    direction    = "y",
    angle        = 0,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    size=3.5, color="black"
  ) +
  geom_label_repel(
    data=dec[cnA==2 & cnB==2 & goodChr==T][order(final_QTL_ID, na.last=F)][PA42qtl==12],
    force   = 10, # do not pull toward data points
    nudge_y      = 75,
    nudge_x       =-5,
    direction    = "y",
    angle        = 0,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    size=3.5, color="black"
  )  +
  geom_label_repel(
    data=dec[cnA==2 & cnB==2 & goodChr==T][order(final_QTL_ID, na.last=F)][PA42qtl%in%c(3,9)],
    force   = 10, # do not pull toward data points
    nudge_y      = 75,
    nudge_x       =-10,
    direction    = "y",
    angle        = 0,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    size=3.5, color="black"
  ) +
  theme(legend.position="none")


  tmrca.plot <-
  ggplot() +
  theme_bw() +
  geom_density(data=tmrca, aes(log10(PMean)), fill="grey", color="grey") +
  geom_point(data=tmrca[!is.na(QTL)], aes(x=log10(PMean), y=0, label=lab), color="black", fill="red", shape=23, size=6) +
  ylab("density") + xlab("log10(TMRCA)") +
  geom_label_repel(
    data=tmrca[!is.na(QTL)],
    aes(x=log10(PMean), y=0, label=lab),
    force_pull   = 0, # do not pull toward data points
    force =100,
    nudge_y      = .4,
    nudge_x       =0,
    direction    = "both",
    angle        = 0,
    hjust        = 0,
    segment.size = 0.2,
    max.iter = 1e4, max.time = 1,
    box.padding=.35, size=3.5
  )


### mega plot

  layout <- "
  AABBCC
  DDDDDD
  EEEEEE
  FFFFFF
  GGGHHH
  GGGHHH
  "

  megaplot <-
  male.plot + epp.plot + f1polar.plot +
  f1.plot +
  poolseq +
  qtlpolar.strip  +
  tmrca.plot + geneexpression.plot +
  plot_layout(design = layout, heights=c(1,1,1,.35,1)) +
  plot_annotation(tag_levels = 'A')

  #ggsave(megaplot, file="~/qtl_mega.pdf")
  ggsave(megaplot, file="~/qtl_mega.png", width=8, height=10)
