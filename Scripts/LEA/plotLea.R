


# scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/Q_samps.Rdata ~/Q_samps.Rdata
# scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/dap.snmf.Rdata ~/dap.snmf.Rdata
# scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/snmf_out.Rdata ~/snmf_out.Rdata
# R

### libraries
  library(LEA)
  library(foreach)
  library(data.table)
  library(ggplot2)
  #library(pophelper)
  library(cowplot); theme_set(theme_cowplot())

### load datda
  #load("~/Q_samps.Rdata")
  #load("~/dap.snmf.Rdata")
  #load("~/snmf_out.v2.Rdata")
  load("~/snmf_out.Rdata")

### which K?
#  ggplot(data=ce.dt, aes(x=k, y=ce, group=run, color=as.factor(run))) + geom_line()
#
#  ce.dt.ag <- ce.dt[,list(mu=median(ce), sd=sd(ce), n=length(ce)), list(k)]
#
#  ggplot() +
#  geom_point(data=ce.dt, aes(x=k, y=ce), fill="grey", alpha=.5) +
#  geom_line(data=ce.dt.ag, aes(x=k, y=mu)) +
#  geom_point(data=ce.dt[k==ce.dt.ag[which.min(mu)]$k][order(ce)][1], aes(x=k, y=ce), color="red")
#
#  ce.dt[which.min(ce)]

### source functions
  source("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/LEA/plotLeaFuncs.R")


### summary plot

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

    plotQ(kr)

  ### hybrids
    kr.ac <- kr[SC%in%c("A", "C")][,c("variable", "SC", "value"), with=F]
    kr.f1 <- kr[AxCF1Hybrid==1][,c("variable", "clone", "value"), with=F]

    #kr.ac.f1 <- merge(kr.ac, kr.f1, by="variable")

    #kr.ac.f1[,exp:=A/2 + C/2]

    ggplot() +
    geom_point(data=kr.f1, aes(x=variable, y=value)) +
    geom_point(data=kr.ac, aes(x=variable, y=value), color="red")



  ### make mega-plot

    q.plot <- plotQ(kr)
    out.plot <- plotQ.sp(kr)


    ce.plot <-   ggplot() +
                 geom_line(data=ce.dt, aes(x=k, y=ce, group=run), color="grey", alpha=.75) +
                 geom_line(data=ce.dt.ag, aes(x=k, y=mu), size=1) +
                 geom_point(data=ce.dt[k==ce.dt.ag[which.min(mu)]$k][order(ce)][1], aes(x=k, y=ce), size=3, color="black") +
                 theme(legend.position="none")

    output.plot <- q.plot / ( out.plot | ce.plot)

    ggsave(output.plot, file="~/LEA_output.new.pdf")












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
