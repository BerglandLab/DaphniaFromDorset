#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(tidyverse)

### files
fn <- c(
      list.files("/project/berglandlab/connor/rohan/out", "bam.summary.txt", full.name=T),
      list.files("/project/berglandlab/connor/rohan/out_all", "bam.summary.txt", full.name=T),
      list.files("/project/berglandlab/connor/rohan/out_other", "bam.summary.txt", full.name=T)
    )


### load data
  rohan <- foreach(fn.i=fn, .combine="rbind")%do%{
    #fn.i<-fn[1]
    tmp <- fread(fn.i, sep="\t", skip=3, header=F)
    tmp.x <- tstrsplit(fn.i, "/")%>%last%>%gsub("_finalmap_mdup.bam.summary.txt", "", .)
    data.table(clone=tmp.x,
               perc_roh=tstrsplit(tmp$V2[4], " ")[[1]]%>%as.numeric,
               len_roh=tstrsplit(tmp$V2[7], " ")[[1]]%>%as.numeric)
  }
  rohan[,clone:=gsub("oo\\.", "", clone)]
  rohan[,clone:=gsub("all.SC.", "", clone)]

### clone metadata file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  sc <- sc[Nonindependent==0 ]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]

### merge
  m <- merge(rohan, sc, by="clone")


### run lengths
  library(doMC)
  registerDoMC(10)
  fn <- c(
        list.files("/project/berglandlab/connor/rohan/out", "mid.hmmrohl", full.name=T),
        list.files("/project/berglandlab/connor/rohan/out_all", "mid.hmmrohl", full.name=T),
        list.files("/project/berglandlab/connor/rohan/out_other", "mid.hmmrohl", full.name=T)
      )

  rl <- foreach(fn.i=fn)%dopar%{
    message(fn.i)
    #fn.i<-fn[1]
    tmp <- fread(fn.i)
    tmp.x <- tstrsplit(fn.i, "/")%>%last%>%gsub("_finalmap_mdup.bam.mid.hmmrohl", "", .)
    tmp[,clone:=tmp.x]
    tmp
  }
  rl <- rbindlist(rl, fill=T)
  rl[,clone:=gsub("oo\\.", "", clone)]
  rl[,clone:=gsub("all.SC.", "", clone)]
  dim(rl)
  rl <- merge(rl, sc, by="clone")
  dim(rl)

  save(m, rl, file="~/rohan_summary.Rdata")


### scp aob2x@rivanna.hpc.virginia.edu:~/rohan_summary.Rdata ~/.

library(ggplot2)
library(data.table)
library(ggrepel)

load("~/rohan_summary.Rdata")

rl.ag <- rl[ROH_LENGTH>=100000, list(nroh=.N, sroh=sum(ROH_LENGTH), mroh=mean(ROH_LENGTH)), list(clone, SC.uniq, population)]
rl.ag.ag <- rl.ag[, list(nroh=mean(nroh), sroh=mean(sroh), mroh=mean(mroh)), list(SC.uniq, population)]


rl.ag.ag[SC.uniq%in%c("A", "C") & population=="D8",lab:=SC.uniq]

rl.ag.ag[is.na(lab), lab:=""]
rl.ag.ag[population=="Dcat", population:="DCat"]

m.ag <- m[,list(perc_roh=mean(perc_roh), len_roh=mean(len_roh)), list(SC.uniq, population)]

t1 <- lm(nroh~sroh, rl.ag.ag[sroh<2.5e6])


summary(lm(sroh~population, rl.ag.ag[population%in%c("D8", "DBunk", "DCat")]))

ggplot(data=rl.ag.ag[population%in%c("D8", "DBunk", "DCat")],
      aes(x=sroh, y=nroh, label=lab, color=as.factor(population), shape=as.factor(population))) +
geom_point() +
geom_label_repel(box.padding=1.5) +
geom_abline(aes(slope=8.181e-06, intercept=2.769e-01))

+
facet_wrap(~population)


ggplot(data=rl.ag, aes(x=sroh, y=nroh, color=as.factor(year), label=SC.uniq)) +
geom_point() +
geom_text_repel() +
facet_wrap(~population)






ggplot(data=m, aes(x=population, y=perc_roh, group=interaction(population, year), fill=as.factor(year))) + geom_boxplot()
ggplot(data=m.ag[population%in%c("D8", "DBunk", "DCat")], aes(x=len_roh, y=perc_roh, color=population, shape=as.factor(population))) + geom_point()
m
