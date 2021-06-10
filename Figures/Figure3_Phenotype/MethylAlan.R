### libraries
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(tidyverse)
        library(lme4)
        library(merTools)
### Load data file
  neo <- fread("FinalNeoNotAdj.csv")
### Alan's version

  neopretrmt <- neo[, c("Clone", "MomJar", "Neonate", "Neo2ndGen", "Isolated", "Clutch1", "Male1", "Female1", "Trmt"), with=FALSE]
  neopretrmtMC <- neopretrmt[Trmt=="MF" | Trmt=="C"][Clone!="SW4" & Clone!="AxoMorse"]
  neopretrmtMC$propmale <- neopretrmtMC$Male1/(neopretrmtMC$Male1+neopretrmtMC$Female1)
  neopretrmtMC$N <- (neopretrmtMC$Male1+neopretrmtMC$Female1)
  neopretrmtMC[,PrePost:=("PreTrmt")]
  neoposttrmt <- neo[, c("Clone", "MomJar", "Neonate", "Neo2ndGen", "Isolated", "Clutch1", "Week3PlusMales", "Week3PlusFemales", "Trmt"), with=FALSE]
  neoposttrmtMC <- neoposttrmt[Trmt=="MF" | Trmt=="C"][Clone!="SW4" & Clone!="AxoMorse"]
  neoposttrmtMC$propmale <- neoposttrmtMC$Week3PlusMales/(neoposttrmtMC$Week3PlusMales+neoposttrmtMC$Week3PlusFemales)
  neoposttrmtMC$N <- (neoposttrmtMC$Week3PlusMales+neoposttrmtMC$Week3PlusFemales)
  neoposttrmtMC[,PrePost:=("PostTrmt")]
  dat <- rbind(neopretrmtMC[,c("Clone", "MomJar", "Trmt", "PrePost", "propmale", "N"), with=F],
              neoposttrmtMC[,c("Clone", "MomJar", "Trmt", "PrePost", "propmale", "N"), with=F])
  dat[grepl("D8515|D8222", Clone), SC:="C"]
  dat[!grepl("D8515|D8222", Clone), SC:="A"]
  #t1 <- glmer(propmale~PrePost*Trmt*SC +(1|Clone:MomJar) , data=dat, family=binomial(), weights=N)
  #t1 <- glmer(propmale~PrePost*Trmt*SC +(1|Clone:MomJar) , data=dat, family=binomial(), weights=N)
  #preds <- predictInterval(t1, n.sims = 1000)
  #pred.dt <- cbind(dat[!is.na(propmale)], preds)
  #pred.dt.ag <- pred.dt[,list(mu.fit=mean(fit), mu.upr=mean(upr), mu.lwr=mean(lwr)), list(Clone, Trmt, PrePost, SC)]
  #pred.dt.ag[,gr:=paste(Clone, Trmt, sep=":")]
  #ggplot(data=pred.dt.ag, aes(x=PrePost, group=gr, y=plogis(mu.fit), color=SC, linetype=Trmt)) +
  #geom_line(position=position_dodge(0.3)) +
  #geom_errorbar(aes(ymin=plogis(mu.lwr), ymax=plogis(mu.upr)), width=.3, position=position_dodge(0.3))

  dat.ag <- dat[,list(propmale=sum(propmale*N, na.rm=T)/sum(N, na.rm=T), N=sum(N, na.rm=T)), list(Clone, Trmt, PrePost, SC)]
  dat.ag[,se:=sqrt((propmale*(1-propmale))/N)]
  dat.ag[,lci:=propmale-1.96*se]
  dat.ag[,uci:=propmale+1.96*se]
  dat.ag[,gr:=paste(Clone, Trmt, sep=":")]
  dat.ag$PrePost <- factor(dat.ag$PrePost, levels=c("PreTrmt", "PostTrmt"))

  ggplot(data=dat.ag, aes(x=PrePost, group=gr, y=propmale, color=SC, linetype=Trmt)) +
    geom_line(position=position_dodge(0.3)) +
    geom_errorbar(aes(ymin=lci, ymax=uci), width=.3, position=position_dodge(0.3)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
