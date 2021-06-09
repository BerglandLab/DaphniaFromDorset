
### libraries
  library(LEA)
  library(foreach)
  library(data.table)
  library(ggplot2)
  library(cowplot); theme_set(theme_cowplot())
  library(patchwork)

### load datda
  load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/SFig_3_hybrid_introgression_plot/snmf_out.Rdata")
  load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/SFig_3_hybrid_introgression_plot/fstat.Rdata")

### process LEA data
  ### source functions
    source("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/LEA/plotLeaFuncs.R")

  ### based on average
    ce.dt.ag <- ce.dt[,list(mu=median(ce), sd=sd(ce), n=length(ce)), list(k)]
    ce.dt.ag[which.min(mu)]
    ce.dt[k==8][which.min(ce)]
    k.i<-8; run.i<-30

  ### absolute lowest
    ce.dt[which.min(ce)]
    k.i<-11; run.i<-21

  ### process data
    kr <- processQ(k=k.i, run=run.i, samp=samps$clone, orderby="set2", n=3)

### Process Dsuite output
  setnames(f, "p-value", "p")
  setnames(f, "f4-ratio", "f4")

  f[data=="orig", pa:=p.adjust(p)]
  f[data=="orig"][P3=="pulicaria"][order(z)][,c('P1', 'P2', 'P3', 'Dstatistic', 'pa', 'f4', 'ABBA', 'BABA'), with=F]
  f[data=="orig"][f4>0.2]

### plot
  # panel A: SNMF individual ancestry plot (structure like plot)
    snmf.ancestry <-
      ggplot(kr, aes(x=ord, y=value, fill=as.factor(variable))) +
      geom_bar(position="stack", stat="identity", width = 1) +
      facet_grid(~set, scale="free_x", space="free") +
      geom_text(data=kr[SC=="A"], aes(x=ord, y=1.1, label="A"), size=5) +
      geom_text(data=kr[SC=="C"], aes(x=ord, y=1.1, label="C"), size=5) +
      theme(strip.text.x = element_text(angle=90),
            panel.spacing = unit(.2, "lines"),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      scale_x_continuous(expand = c(0,0)) +
      xlab("") + ylab("") +
      guides(fill=FALSE)

  ### panels B & C: Coefficients for Dpuli & Dobt ancestry
    kr[SC=="A",col:="A"]
    kr[SC=="C",col:="C"]
    setorderv(kr, "col")

    puli.k <- kr[Species=="pulicaria"][which.max(value)]$variable
    obt.k <- kr[Species=="obtusa"][which.max(value)]$variable

    puli.plot <-
      ggplot(data=kr[variable==puli.k][!is.na(set)][grepl("pulex", set)], aes(x=set, y=log10(value), color=col)) +
      geom_jitter(width = 0.25) +
      coord_flip() +
      ylab("log10(% ancestry):\nD. pulicaria") + xlab("") +
      scale_x_discrete(labels= c("DRamp", "DOil", "DCat", "D8", "DBunk", "D10", "W1", "W2")) +
      theme(legend.position="none")

    obt.plot <-
      ggplot(data=kr[variable==obt.k][!is.na(set)][grepl("pulex", set)], aes(x=set, y=log10(value), color=col)) +
      geom_jitter(width = 0.25) +
      coord_flip() +
      ylab("log10(% ancestry):\nD.obtusa") + xlab("") +
      scale_x_discrete(labels= NULL) +
      theme(legend.position="none")

  ### panel D: optimal K
    ce.plot <-
      ggplot() +
      geom_line(data=ce.dt, aes(x=k, y=ce, group=run), color="grey", alpha=.75) +
      geom_line(data=ce.dt.ag, aes(x=k, y=mu), size=1) +
      geom_point(data=ce.dt[k==ce.dt.ag[which.min(mu)]$k][order(ce)][1], aes(x=k, y=ce), size=3, color="black") +
      theme(legend.position="none") +
      ylab("Cross Entropy")


  ### panel D: Dsutie analysis
    f4 <-
      ggplot() +
      geom_point(data=f[data=="orig"][P3=="pulicaria"],
                  aes(x=Dstatistic, y=log10(f4),
                     shape=as.factor(pa<.05)), size=4) +
      theme(legend.position = c(0.8, 0.2)) +
      ylab("log10(f4-ratio)") +
      scale_shape_manual(values = c(1,13), name="", labels=c("pa > 0.05", "pa < 0.05"))


### merge
layout <- "
AAAAAA
BCCDDD
##EEEE
"

big.plot <-
snmf.ancestry + ce.plot + puli.plot + obt.plot + f4 +
plot_annotation(tag_levels = 'A', theme=theme(plot.tag = element_text(size = 19))) +
plot_layout(design = layout)

### output
  ggsave(big.plot, file="/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/SFig_3_hybrid_introgression_plot/LEA_output.new.pdf", w=7, h=9)












layout <- "
AAAAAAAA
##BBBBCC
"




  +
  plot_layout(design = layout)


plot.q1 <- plotQ(processQ(k=1, samp=samps$clone, orderby="set"))
plot.q2 <- plotQ(processQ(k=2, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q3 <- plotQ(processQ(k=3, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q4 <- plotQ(processQ(k=4, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q5 <- plotQ(processQ(k=5, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q6 <- plotQ(processQ(k=6, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())
plot.q7 <- plotQ(processQ(k=7, samp=samps$clone, orderby="set")) + theme(strip.background = element_blank(),strip.text.x = element_blank())

plot.q1 / plot.q2 / plot.q3 / plot.q4 / plot.q5 / plot.q6 / plot.q7



hist(k[variable.x=="V6"]$value.x)
k[variable.x=="V6"][value.x>.2]


hist(log10(k[variable.x=="V3"]$value.x), breaks=1000)
hist(log10(k[variable.x=="V6"]$value.x), breaks=1000)
