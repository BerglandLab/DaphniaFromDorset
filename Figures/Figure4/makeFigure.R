### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(patchwork)
  library(lme4)


### setwd
  setwd("/Users/kbkubow/Documents/GitHub")

### load F1 cross data (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")

### load in PoolSeq data (made by `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLSeqR.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")

### load in F1 mapping data
  #load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")
  #load("~/lme4qtl_output.AxC.long.Rdata")
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.SUMMARIZED.Rdata")

### load in overlap tests
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/overlap_perm.Rdata")

### Load in PA42 chr key
  PA42 <- fread("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/5kbchrassignHiCnew.csv")

### Load in qtl_polar
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/qtl_polar.Rdata")

### male proportion plot
  male$gr <- ifelse(male$SC=="selfedA", "A", ifelse(male$SC=="B", "B", male$gr))
  male.ag <- male[!is.na(gr),list(propmale=sum(Males)/sum(NewTotal),
                                  N=sum(NewTotal)),
                    list(clone, gr)]
  male.ag[,se:=sqrt((propmale*(1-propmale))/N)]

  male.ag[,lci:=propmale-1.96*se]
  male.ag[,uci:=propmale+1.96*se]

  male.ag$gr <- factor(male.ag$gr, levels=c("A", "AxC", "C", "CxC", "B"))

  male.plot <- ggplot(data=male.ag, aes(x=gr, y=propmale, color=gr)) +
  geom_linerange(aes(ymin=lci, ymax=uci), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Male proportion") +theme_bw()

#varmaleA <- var(male.ag$propmale[male.ag$gr=="A"])
#varmaleC <- var(male.ag$propmale[male.ag$gr=="C"])
#sampvarmaleA <- varmaleA/5
#sampvarmaleC <- varmaleC/2
#varmaleAxC <- var(male.ag$propmale[male.ag$gr=="AxC"])
#meanA <- mean(male.ag$propmale[male.ag$gr=="A"])
#meanC <- mean(male.ag$propmale[male.ag$gr=="C"])
#nummalelociAxC <- ((meanA-meanC)*(meanA-meanC)-sampvarmaleA-sampvarmaleC)/(8*varmaleAxC)

#varmaleCxC <- var(male.ag$propmale[male.ag$gr=="CxC"])
#nummalelociCxC<- ((male.ag$propmale[male.ag$clone=="April_2017_D8_222"]-
#    male.ag$propmale[male.ag$clone=="May_2017_D8_515"])*(male.ag$propmale[male.ag$clone=="April_2017_D8_222"]-
#        male.ag$propmale[male.ag$clone=="May_2017_D8_515"])-envvarmaleC)/(8*varmaleCxC)

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

  epp.ag[,seemb:=sdEMB/(sqrt(sampsize))]
  epp.ag[,lciemb:=meanEMB-1.96*seemb]
  epp.ag[,uciemb:=meanEMB+1.96*seemb]

  epp.ag[,seepp:=sdEpp/(sqrt(sampsize))]
  epp.ag[,lciepp:=meanEpp-1.96*seepp]
  epp.ag[,uciepp:=meanEpp+1.96*seepp]


  epp.ag$gr <- factor(epp.ag$gr, levels=c("A", "AxC", "C", "CxC", "B"))

  epp.plot <- ggplot(data=epp.ag, aes(x=gr, y=fillrate, color=gr)) +
  geom_linerange(aes(ymin=lci, ymax=uci), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Ephippial fill rate") +
  theme_bw()

  emb.plot <- ggplot(data=epp.ag, aes(x=gr, y=meanEMB, color=gr)) +
  geom_linerange(aes(ymin=lciemb, ymax=uciemb), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Embryo production") +
  theme_bw()

  eppnum.plot <- ggplot(data=epp.ag, aes(x=gr, y=meanEpp, color=gr)) +
  geom_linerange(aes(ymin=lciepp, ymax=uciepp), position = position_jitter(seed = 123, width =0.2), size=.25, alpha=.75, color="black") +
  geom_point(position = position_jitter(seed = 123, width =0.2), size=2) +
  theme(legend.position="none") +
  xlab("Clone or Cross Type") + ylab("Ephippia production") +
  theme_bw()

  varmaleA <- var(epp.ag$fillrate[epp.ag$gr=="A"])
  varmaleC <- var(epp.ag$fillrate[epp.ag$gr=="C"])
  sampvarmaleA <- varmaleA/5
  sampvarmaleC <- varmaleC/2
  varmaleAxC <- var(epp.ag$fillrate[epp.ag$gr=="AxC"])
  meanA <- mean(epp.ag$fillrate[epp.ag$gr=="A"])
  meanC <- mean(epp.ag$fillrate[epp.ag$gr=="C"])
  nummalelociAxC <- ((meanA-meanC)*(meanA-meanC)-sampvarmaleA-sampvarmaleC)/(8*varmaleAxC)


### trait correlation plot
  setkey(epp.ag, clone, gr)
  setkey(male.ag, clone, gr)
  f1 <- merge(epp.ag, male.ag)

  propmalefillrate <- ggplot(data=f1, aes(x=propmale, y=fillrate, color=gr)) + geom_point() +
    theme_bw() + xlab("Proportion male") + ylab("Embryos per ephippia")
  eppemb <- ggplot(data=f1, aes(x=meanEpp, y=meanEMB, color=gr)) + geom_point() +
    geom_abline(intercept=0, slope=2, linetype=2) + theme_bw() + xlab("Mean ephippia") +
    ylab("Mean embryos")

  propmaleepp <- ggplot(data=f1, aes(x=propmale, y=meanEpp, color=gr)) + geom_point() + theme_bw()
  propmaleemb <- ggplot(data=f1, aes(x=propmale, y=meanEMB, color=gr)) + geom_point() + theme_bw()
  eppfillrate <- ggplot(data=f1, aes(x=meanEpp, y=fillrate, color=gr)) + geom_point() + theme_bw()


  fillrateemb <- ggplot(data=f1, aes(x=fillrate, y=meanEMB, color=gr)) + geom_point() + theme_bw()

  F1phenoplot <- propmalefillrate + propmaleepp + propmaleemb + eppfillrate + eppemb +
    fillrateemb + plot_layout(nrow=2, byrow=TRUE, guides='collect') + plot_annotation(tag_levels = 'A')

  cor.test(f1$propmale, f1$fillrate)

    #Pearson's product-moment correlation

    #data:  f1$propmale and f1$fillrate
    #t = 2.607, df = 56, p-value = 0.01168
    3alternative hypothesis: true correlation is not equal to 0
    395 percent confidence interval:
    #0.07725061 0.54128251
    #sample estimates:
     #cor
     #0.328982

  cor.test(f1$meanEpp, f1$fillrate)

   	#Pearson's product-moment correlation

    #data:  f1$meanEpp and f1$fillrate
    #t = 2.1741, df = 56, p-value = 0.03394
    #alternative hypothesis: true correlation is not equal to 0
    #95 percent confidence interval:
    #0.02230105 0.50117028
    #sample estimates:
    #     cor
    #0.2789898

cor.test(f1$meanEpp, f1$meanEMB)

    #Pearson's product-moment correlation

    #data:  f1$meanEpp and f1$meanEMB
    #t = 17.791, df = 58, p-value < 2.2e-16
    #alternative hypothesis: true correlation is not equal to 0
    #95 percent confidence interval:
    # 0.8680078 0.9511919
    #sample estimates:
    #     cor
    #0.9193089

cor.test(f1$fillrate, f1$meanEMB)

   	#Pearson's product-moment correlation

    #data:  f1$fillrate and f1$meanEMB
    #t = 4.5838, df = 56, p-value = 2.604e-05
    #alternative hypothesis: true correlation is not equal to 0
    #95 percent confidence interval:
    #0.3052126 0.6878274
    #sample estimates:
    #     cor
    #0.5223304


  # Correlation between propmale and meanEpp not significant
  # Correlation between propmale and meanEMB also not significant


  t1 <- lm(meanEMB ~ meanEpp, data=f1)
  t2 <- lm(meanEMB ~ meanEpp + propmale, data=f1)

  t1 <- lm(meanEMB ~ meanEpp, data=f1[gr=="AxC"])
  t2 <- lm(meanEMB ~ meanEpp + propmale, data=f1[gr=="AxC"])

  t1 <- lm(meanEMB ~ meanEpp, data=f1[gr=="CxC"])
  t2 <- lm(meanEMB ~ meanEpp + propmale, data=f1[gr=="CxC"])

  anova(t1, t2)


  t3 <- lm(fillrate ~ meanEpp, data=f1)
  t4 <- lm(fillrate ~ meanEpp + propmale, data=f1)

  t3 <- lm(fillrate ~ meanEpp, data=f1[gr=="AxC"])
  t4 <- lm(fillrate ~ meanEpp + propmale, data=f1[gr=="AxC"])

  t3 <- lm(fillrate ~ meanEpp, data=f1[gr=="CxC"])
  t4 <- lm(fillrate ~ meanEpp + propmale, data=f1[gr=="CxC"])

  anova(t3, t4)



  t1emb <- lmer(EmbDividebyTot~1+(1|Week)+SC/SCrep/Replicate, data=meso7AC, weights=TotalInd)
  t2emb <- lmer(EmbDividebyTot~1+(1|Week), data=meso7AC, weights=TotalInd)


### Pool seq results.
    setnames(peaks, "CHROM", "chr")
    setnames(gprime, "CHROM", "chr")
    #setnames(peaks.groups, "CHROM", "chr")
    #setnames(gprime.groups, "CHROM", "chr")

    PA42$hybridchr <- paste("Chr", PA42$PA42chr, "\n", paste("Scaff", "\n", tstrsplit(PA42$chr, "_")[[2]], sep=""), sep="")
    PA42sub <- PA42[, c("chr", "PA42chr", "hybridchr")]
    setkey(PA42sub, chr)
    setkey(peaks, chr)
    setkey(gprime, chr)
    mpeaks <- merge(peaks, PA42sub)
    mgprime <- merge(gprime, PA42sub)

    #peaks <- rbind(peaks, peaks.groups, fill=T)
    #gprime <- rbind(gprime, gprime.groups, fill=T)

    #peaks[,rep:=gsub("NA", "", paste(rep, group, sep=""))]
    #gprime[,rep:=gsub("NA", "", paste(rep, group, sep=""))]

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
    poolseq <- ggplot() +
    #geom_vline(data=mpeaks, aes(xintercept=posMaxGprime), color="black") +
    geom_line(data=mgprime, aes(x=pos, y=Gprime.y, group=rep, color=as.factor(rep)), size=1.5) +
    geom_segment(data=mpeaks, aes(x=posMaxGprime, xend=posMaxGprime,
                              y=maxGprime.y+0.15*ifelse(rep==1, 1, -1), yend=5*ifelse(rep==1, 1, -1))) +
    geom_hline(data=mgprime[,list(minG=min(Gprime[qvalue<=0.05])), list(rep)],
           aes(yintercept=minG*ifelse(rep==1, 1, -1))) +
    geom_text(data=mpeaks, aes(x=posMaxGprime, y=5.5*ifelse(rep==1, 1, -1), label=c(1:14)), size=3) +
    facet_grid(~hybridchr, scales="free_x", space = "free_x") +
    theme(legend.position = "none") +
    theme_bw() +
    scale_x_continuous(breaks=seq(from=1, to=13, by=2)*10^6, labels=seq(from=1, to=13, by=2)) +
    theme(strip.text.x = element_text(size = 8)) + ylab("Gprime") + xlab("Position (Mb)") +
    labs(color="Replicate")

    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +

### F1 mapping results
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

  f1.plot <- ggplot() +
  geom_vline(data=mpeaks, aes(xintercept=posMaxGprime, linetype=as.factor(rep))) +
  geom_line(data=mo.ag.plot, aes(x=pos, y=-log10(p.aov), color=hybridchr), size=1.5) +
  geom_point(data=mo.ag.plot[p.aov<=p.aov.thr], aes(x=pos, y=-log10(p.aov)), size=.5) +
  geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.aov.thr), linetype=as.factor(q))) +
  facet_grid(term~hybridchr, scales="free_x", space = "free_x") +
  theme_bw() + theme(strip.text.x = element_text(size = 8)) +
  scale_x_continuous(breaks=seq(from=1, to=13, by=2)*10^6, labels=seq(from=1, to=13, by=2)) +
  theme(legend.position = "none") + xlab("Position (Mb)")

### Qtl polar plot

  qtlcoding <- data.table(qtl=c(1:14), PA42qtl=c(3, 4, 2, 1, 9, 10, 11, 13, 14, 12, 8, 6, 7, 5))

  mqtl.polar <- merge(qtl.polar, qtlcoding, by="qtl")
  qtl.polar <- mqtl.polar
  setkey(mqtl.polar, PA42qtl)

  qtl.polar$geno <- ifelse(qtl.polar$n.male_male > qtl.polar$n.male_pe &
    qtl.polar$n.male_male > qtl.polar$n.pe_pe, "male_male",
    ifelse(qtl.polar$n.male_pe > qtl.polar$n.male_male &
    qtl.polar$n.male_pe > qtl.polar$n.pe_pe, "male_pe",
    ifelse(qtl.polar$n.pe_pe > qtl.polar$n.male_male &
    qtl.polar$n.pe_pe > qtl.polar$n.male_pe,"pe_pe", "other")))

  qtlsub <- qtl.polar[SC.uniq=="C" & year=="2017" & pond.y=="D8" | SC.uniq=="A" &
    year=="2017" & pond.y=="D8"| SC.uniq=="B" & year=="2018"]

  qtlcounts <- qtlsub[, .N, by=list(SC.uniq, geno)]
  Cpepe <- data.table(SC.uniq="C", geno="pe_pe", N=0)
  qtlcountsB <- rbind(qtlcounts, Cpepe)

  qtlpolar <- ggplot(data=qtlcountsB, aes(fill=geno, x=SC.uniq, y=N)) + geom_bar(position="stack",
    stat="identity") + theme_bw() + xlab("Clonal lineage") + labs(fill="Genotype")

  ggplot(data=qtlsub, aes(x=SC.uniq, y=PA42qtl, fill=geno)) + geom_tile() +
    xlab("Clonal lineage") + ylab("QTL") + labs(fill="Genotype")

  qtlsub$qtl <- factor(qtlsub$qtl, levels=c(2, 5, 13, 1, 9, 11, 14, 4, 12, 10, 7, 3, 6, 8))

  qtlcounts$geno <- factor(qtlcounts$geno, levels=c("male_male", "male_pe", "pe_pe"))

  qtlpolar <- ggplot(data=qtlcounts, aes(x=geno, y=N, fill=geno)) + geom_bar(stat="identity") +
    facet_grid(~SC.uniq) + theme_bw() + theme(legend.position = "none") + xlab("Genotype")

### overlap plot
  overlap.plot <- ggplot() +
  geom_histogram(data=overlap.perm[perm!=0], aes(z)) +
  geom_vline(data=overlap.perm[perm==0], aes(xintercept=z)) +
  facet_grid(cross~pheno)

### F1 qtl polar plot

### load
  library(data.table)
  library(ggplot2)

  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/f1_pool_polar.Rdata")
  load("f1_pool_polar.Rdata")

  f1.pool.mergeB <- merge(f1.pool.merge, qtlcoding, by="qtl")
  f1.pool.merge <- f1.pool.mergeB

  m.ag <- f1.pool.merge[,list(propmale=mean(propmale), sd=sd(propmale)), list(qtl, geno, PA42qtl)]
  sampsize <- length(unique(f1.pool.merge$clone))

  m.ag[,se:=sd/(sqrt(sampsize))]
  m.ag[,lci:=propmale-1.96*se]
  m.ag[,uci:=propmale+1.96*se]

  m2.ag <- m.ag[,list(propmale=mean(propmale), sd=sd(propmale)), list(geno)]

  ggplot() +
  geom_point(data=f1.pool.merge, aes(x=geno, y=propmale, color=gr), size=.75, alpha=.95) +
  geom_point(data=m.ag, aes(x=geno, y=propmale), size=2, color="red") +
  facet_wrap(~qtl)

  ggplot(data=m.ag[qtl!=10], aes(x=geno, y=propmale, group=qtl, linetype=as.factor(qtl))) + geom_line() +
    geom_line(data=m.ag[qtl=="10"], aes(x=geno, y=propmale, group=qtl, color="QTL_12"), size=2) +
    theme_bw() + theme(legend.position = "none")

    f1.pool.mergesub <- f1.pool.merge[, c("clone", "qtl", "geno"), with=TRUE]
    malesqtl <- merge(f1.pool.mergesub, male, by="clone", allow.cartesian=TRUE)

    anovas <- foreach(q=1:14, .combine="rbind")%do%{
      f1sub <- malesqtl[qtl==q]
      f1subout_a <- glmer(propmale~(1|clone) + (1|Replicate), data=f1sub, family=binomial(), weights=NewTotal)
      f1subout_b <- glmer(propmale~geno+(1|clone) + (1|Replicate), data=f1sub, family=binomial(), weights=NewTotal)
      aovcompare <- anova(f1subout_a, f1subout_b)
      p <- aovcompare[[8]][2]
      tmp <- data.table(qtl=q, p=p)
      tmp
      }

    m.agsig <- merge(m.ag, anovas, by="qtl")

    F1polar <- ggplot(data=m.agsig[p <= 0.05], aes(x=geno, y=propmale, group=qtl)) + geom_line() +
      theme_bw() + theme(legend.position = "none") + xlab("Genotype") + ylab("Proportion male") +
      geom_errorbar(aes(ymin=propmale-2*se, ymax=propmale+2*se), width=0.2) +
      facet_wrap(~PA42qtl)

    F1polarqtl10 <- ggplot(data=m.agsig[qtl==10], aes(x=geno, y=propmale, group=qtl, linetype=as.factor(qtl))) + geom_line() +
      geom_errorbar(aes(ymin=propmale-2*se, ymax=propmale+2*se), width=0.2) +
      theme_bw() + theme(legend.position = "none") + xlab("Genotype") + ylab("Proportion male") +
      geom_linerange(aes(ymin=lci, ymax=uci))

    F1polarqtl5and12 <- ggplot(data=m.agsig[PA42qtl==5 | PA42qtl==12], aes(x=geno, y=propmale, group=qtl, linetype=as.factor(PA42qtl))) + geom_line() +
      geom_errorbar(aes(ymin=propmale-2*se, ymax=propmale+2*se), width=0.2) +
      theme_bw() + xlab("Genotype") + ylab("Proportion male") + labs(linetype="QTL")


    #### abundance genotype
      ab <- qtl.polar[,list(geno=c("male_male", "male_pe", "pe_pe")[which.max(c(n.male_male, n.male_pe, n.pe_pe))],
                            size=N[1]), list(clone, pond=pond.y, qtl, PA42qtl, sc=SC.uniq, year=year)]

      abr <- ab[,list(geno=rep(geno, size)), list(clone, pond, sc, qtl, PA42qtl, year)]

      abrf <- abr[,list(male_freq=(2*sum(geno=="male_male") + sum(geno=="male_pe"))/(2*length(geno)), n=2*length(geno)), list(pond, year, qtl, PA42qtl)]
      abrf[,se:=male_freq*(1-male_freq)/sqrt(n)]
      #abrf[,pond:=factor(pond, levels=c("DBUNK", "D8", "DCAT"))]

      good <- ggplot(data=abrf[pond%in%c("D8", "DBUNK", "DCAT")][year>2016], aes(x=as.factor(year), y=male_freq, group=pond, color=pond)) +
      geom_errorbar(aes(ymin=male_freq-2*se, ymax=male_freq+2*se), width=0.5) +
      geom_point() + geom_line() + facet_grid(~PA42qtl) +
      xlab("Year") + ylab("Fequency of male allele") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      goodqtl10 <- ggplot(data=abrf[pond%in%c("D8", "DBUNK", "DCAT")][year>2016 & PA42qtl==12], aes(x=as.factor(year), y=male_freq, group=pond, color=pond)) +
      geom_errorbar(aes(ymin=male_freq-2*se, ymax=male_freq+2*se), width=0.2) +
      geom_point() + geom_line() +
      xlab("Year") + ylab("Freq male allele") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      goodqtl10 <- ggplot(data=abrf[pond%in%c("D8", "DBUNK", "DCAT")][year>2016 & PA42qtl==12 |
      year>2016 & PA42qtl==5], aes(x=as.factor(year), y=male_freq, group=pond, color=pond)) +
      geom_errorbar(aes(ymin=male_freq-2*se, ymax=male_freq+2*se), width=0.2) +
      geom_point() + geom_line() + facet_wrap(~PA42qtl)
      xlab("Year") + ylab("Freq male allele") + theme_bw()


### mega plot

phenoplot <- propmalefillrate + eppemb +
  plot_layout(nrow=1, byrow=TRUE, guides='collect') + plot_annotation(tag_levels = 'A')

mapping <- poolseq + f1.plot +
  plot_layout(nrow=2, byrow=TRUE)

confirm <- qtlpolar + F1polarqtl10 + goodqtl10 +
  plot_layout(nrow=1, byrow=TRUE)

totalfig4 <- poolseq /phenoplot / f1.plot / confirm + plot_layout(heights = c(2, 2.5, 2, 2)) +
  plot_annotation(tag_levels = 'A')




eppnum.plot

layout <- "
ABC
DDD
EEE
"

bigplot <- male.plot + epp.plot + emb.plot + poolseq + f1.plot +
plot_annotation(tag_levels = 'A', theme=theme(plot.tag = element_text(size = 19))) +
plot_layout(design = layout)

ggsave(bigplot, file="~/big_qtl.png", h=8, w=11)














###
