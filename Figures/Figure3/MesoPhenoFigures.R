#!/usr/bin/env Rscript

### libraries
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(tidyverse)
        library(lme4)
        library(cowplot)
        library(merTools)
        library(patchwork)

### Load data file

  meso <- fread("SampleCountsMeso2019Week7.csv")
  meso$CloneRep <- paste(meso$Clone, meso$Replicate, sep="_")
  meso$SC <- ifelse(meso$Clone=="D8179" | meso$Clone=="D8349", "A", ifelse(
    meso$Clone=="D8222" | meso$Clone=="D8515", "C", "AxC"
  ))

### Cleaning up and adding new variables
  meso$Clone <- as.factor(meso$Clone)
  meso$TotalPop <- meso$TotalInd*15
  meso[is.na(JuvMale),JuvMale:=0]
  meso[is.na(AdMale),AdMale:=0]
  meso[is.na(EstimatedEmbryosTot),EstimatedEmbryosTot:=0]
  meso$totmale <- meso$JuvMale + meso$AdMale
  meso$EmbDividebyTot <- meso$EstimatedEmbryosTot/meso$TotalPop
  meso$EppDividebyTot <- meso$LooseEppTotal/meso$TotalPop
  meso$SCrep <- ifelse(meso$SC=="A" & meso$Clone=="D8179", "1", ifelse(meso$SC=="C" & meso$Clone=="D8222", "1", "2"))
  meso[is.na(MomPE_5),MomPE_5:=0]
  meso[is.na(MomPE_10),MomPE_10:=0]
  meso[is.na(MomPE_15),MomPE_15:=0]
  meso[is.na(MomPE_20),MomPE_20:=0]
  meso[is.na(MomwEpp),MomwEpp:=0]
  meso[is.na(Barren),Barren:=0]
  meso$totmomPE <- meso$MomPE_5 + meso$MomPE_10 + meso$MomPE_15 + meso$MomPE_20
  meso$totalADfemale <- meso$totmomPE + meso$MomwEpp + meso$Barren
  meso$propPEtotalind <- meso$totmomPE/meso$TotalInd
  meso$propPEmoms <- meso$totmomPE/meso$totalADfemale
  meso$propMomEpp <- meso$MomwEpp/meso$TotalInd
  meso$propfill <- meso$EstimatedEmbryosTot/(meso$LooseEppTotal*2)

### Restrict to weeks 1-7 and SCs A and C (not looking at AxC)

  meso7AC <- meso[Week<8 & SC=="A" | Week<8 & SC=="C"]

### Figures
  week7ACtotalind <- ggplot(data=meso7AC, aes(x=Week, y=TotalPop, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("# Individuals")


  week7ACmomPEnormtotal <- ggplot(data=meso7AC, aes(x=TotalInd, y=propPEtotalind, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Asexual Females") + xlab("Total Individuals")

  week7ACmomPEnormCumInd <- ggplot(data=meso7AC, aes(x=CumInd, y=propPEtotalind, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Asexual Females") + xlab("Cumulative Individuals")

  week7ACmomPEnormWeek <- ggplot(data=meso7AC, aes(x=Week, y=propPEtotalind, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Asexual Females")


  week7ACLooseEppnormtotal <- ggplot(data=meso7AC, aes(x=TotalInd, y=EppDividebyTot, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Loose Epphippia Per Individual") + xlab("Total Individuals")

  week7ACLoseeEppnormCumInd <- ggplot(data=meso7AC, aes(x=CumInd, y=EppDividebyTot, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Loose Epphippia Per Individual") + xlab("Cumulative Individuals")

  week7ACLoseeEppnormWeek <- ggplot(data=meso7AC, aes(x=Week, y=EppDividebyTot, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Loose Epphippia/Total Individuals")


  week7ACPropMalenormtotal <- ggplot(data=meso7AC, aes(x=TotalInd, y=PropTotMale, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Male") + xlab("Total Individuals")

  week7ACPropMalenormCumInd <- ggplot(data=meso7AC, aes(x=CumInd, y=PropTotMale, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Male") + xlab("Cumulative Individuals")

  week7ACPropMalenormWeek <- ggplot(data=meso7AC, aes(x=Week, y=PropTotMale, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Male")


  week7ACEmbnormtotal <- ggplot(data=meso7AC, aes(x=TotalInd, y=EmbDividebyTot, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Embryos/Total Individuals") + xlab("Total Individuals")

  week7ACEmbnormCumInd <- ggplot(data=meso7AC, aes(x=CumInd, y=EmbDividebyTot, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Embryos/Total Individuals") + xlab("Cumulative Individuals")

  week7ACEmbnormWeek <- ggplot(data=meso7AC, aes(x=Week, y=EmbDividebyTot, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Embryos/Total Individuals")

  week7ACEmbWeek <- ggplot(data=meso7AC, aes(x=Week, y=EstimatedEmbryosTot, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("# Sexual Embryos")

  from <- 250
  to <- 700

  week7ACEmbWeekB <- ggplot(data=meso7AC, aes(x=Week, y=log10(EstimatedEmbryosTot), color=SC,
    group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("log10(# Sexual Embryos)") +
    scale_y_continuous(breaks=c(0, 1, 2, 3), labels=c("0", "10", "100", "1000"))



  week7ACMomwEppnormtotal <- ggplot(data=meso7AC, aes(x=TotalInd, y=propMomEpp, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Sexual Females") + xlab("Total Individuals")

  week7ACMomwEppnormCumInd <- ggplot(data=meso7AC, aes(x=CumInd, y=propMomEpp, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Sexual Females") + xlab("Cumulative Individuals")

  week7ACMomwEppnormWeek <- ggplot(data=meso7AC, aes(x=Week, y=propMomEpp, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Sexual Females")


  week7ACMomwPropFilltotal <- ggplot(data=meso7AC, aes(x=TotalInd, y=propfill, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Sexual Females") + xlab("Total Individuals")

  week7ACMomwPropFillCumInd <- ggplot(data=meso7AC, aes(x=CumInd, y=propfill, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Sexual Females") + xlab("Cumulative Individuals")

  week7ACMomwPropFillWeek <- ggplot(data=meso7AC, aes(x=Week, y=propfill, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Sexual Females")


  week7ACEstEmb <- ggplot(data=meso7AC, aes(x=Week, y=EmbDividebyTot, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACPropTotMale <- ggplot(data=meso7AC, aes(x=Week, y=PropTotMale, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

### Arranging in grid layout

  totalCumInd <- week7ACtotalind + week7ACmomPEnormWeek + week7ACLoseeEppnormWeek + week7ACPropMalenormWeek +
        week7ACEmbWeek + week7ACmomPEnormCumInd + week7ACLoseeEppnormCumInd + week7ACPropMalenormCumInd +
        plot_layout(nrow = 2, byrow = TRUE, guides = 'collect')

  totalCumIndmomwEpp <- week7ACtotalind + week7ACmomPEnormWeek + week7ACMomwEppnormWeek + week7ACPropMalenormWeek +
        week7ACEmbWeek + week7ACmomPEnormCumInd + week7ACMomwEppnormCumInd + week7ACPropMalenormCumInd +
        plot_layout(nrow = 2, byrow = TRUE, guides = 'collect')

### Add in methylfarnesoate figure
### Load data file
  neo <- fread("/Users/kbkubow/Box Sync/Daphnia/Mesocosms2019/MethylFanesoate/FinalNeoNotAdj.csv")
### Alan's version

  neopretrmt <- neo[, c("Clone", "MomJar", "Neonate", "Neo2ndGen", "Isolated", "Clutch1", "Male1", "Female1", "Trmt"), with=FALSE]
  neopretrmtMC <- neopretrmt[Trmt=="MF" | Trmt=="C"][Clone!="SW4" & Clone!="AxoMorse"]
  neopretrmtMC$propmale <- neopretrmtMC$Male1/(neopretrmtMC$Male1+neopretrmtMC$Female1)
  neopretrmtMC$N <- (neopretrmtMC$Male1+neopretrmtMC$Female1)
  neopretrmtMC[,PrePost:=("Pre")]
  neoposttrmt <- neo[, c("Clone", "MomJar", "Neonate", "Neo2ndGen", "Isolated", "Clutch1", "Week3PlusMales", "Week3PlusFemales", "Trmt"), with=FALSE]
  neoposttrmtMC <- neoposttrmt[Trmt=="MF" | Trmt=="C"][Clone!="SW4" & Clone!="AxoMorse"]
  neoposttrmtMC$propmale <- neoposttrmtMC$Week3PlusMales/(neoposttrmtMC$Week3PlusMales+neoposttrmtMC$Week3PlusFemales)
  neoposttrmtMC$N <- (neoposttrmtMC$Week3PlusMales+neoposttrmtMC$Week3PlusFemales)
  neoposttrmtMC[,PrePost:=("Post")]
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
  dat.ag$PrePost <- factor(dat.ag$PrePost, levels=c("Pre", "Post"))

  Methyl <- ggplot(data=dat.ag, aes(x=PrePost, group=gr, y=propmale, color=SC, linetype=Trmt)) +
    geom_line(position=position_dodge(0.3)) +
    geom_errorbar(aes(ymin=lci, ymax=uci), width=.3, position=position_dodge(0.3)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    xlab("Pre/Post Treatment") + ylab("Prop Male")


### OK, let's try to get in the 1 liter and 250 ml data
### Load in data files
    onelitermales <- fread("/Users/kbkubow/Box Sync/Daphnia/A_B_1liters/MaleCensusAC.csv")
    grace250spring <- fread("/Users/kbkubow/Box Sync/Daphnia/A_B_1liters/Grace250Sp2020.csv")
    grace250spring$SC <- ifelse(grace250spring$Clone=="D8179" | grace250spring$Clone=="D8349", "A", "C")
    grace250spring$momsplusslide <- grace250spring$Males_4 + grace250spring$Moms_4 + grace250spring$FemalesSlide_4
    grace250spring$propmale <- grace250spring$Males_4/grace250spring$momsplusslide
    grace250spring$maleyesno <- ifelse(grace250spring$Males_4 > 0, 1, 0)
    grace250spring$count <- c(1)

### Calculate proportion males
    onelitermales$propmale <- onelitermales$Males/onelitermales$NewTotal

    grace250spring.ag <- grace250spring[,list(propmale=sum(maleyesno, na.rm=T)/sum(count, na.rm=T)), list(Clone, SC)]



### Combine datasets
  #Let's take just the endpoints of the mesocosm and 250 ml data
      mesoWeek7 <- meso7AC[Week==7]
  # Now simplify
      mesoWeek7sub <- mesoWeek7[, c("SC", "Clone", "Replicate", "PropTotMale"), with=TRUE]
      colnames(mesoWeek7sub) <- c("SC", "Clone", "Replicate", "propmale")
      mesoWeek7sub$exp <- c("Mesocosm")
  # Now simplify one liters
      onelitermalessub <- onelitermales[, c("SC", "Clone", "Replicate", "propmale"), with=TRUE]
      onelitermalessub$exp <- c("OneLiters")
### Add in Grace's data here
  # Now to combine
      totalmale <- rbind(mesoWeek7sub, onelitermalessub)

### Now to plot
      totalmale[,Clone := factor(Clone, levels=c("D8179", "D8349", "D8222", "D8515"))]
      MesoandOneLiter <- ggplot(data=totalmale, aes(x=Clone, y=propmale, color=SC)) + geom_point() + facet_wrap(~exp) + theme_cowplot(12) + ylim(0,0.45)
      Meso <- ggplot(data=totalmale[exp=="Mesocosm"], aes(x=Clone, y=propmale, color=SC)) + geom_point() +
        ylim(0,0.45) + ylab("Prop Male") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        theme(legend.position = "none")

      OneLiter <- ggplot(data=totalmale[exp=="OneLiters"], aes(x=Clone, y=propmale, color=SC)) + geom_point() +
        ylim(0,0.45) + ylab("Prop Male") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        theme(legend.position = "none")


### Now let's add in Graces' 250 ml - box plot may be better here
      grace250spring[,Clone := factor(Clone, levels=c("D8179", "D8349", "D8222", "D8515"))]
      Grace250 <- ggplot(data=grace250spring, aes(x=Clone, y=propmale, color=SC, group=Clone)) +
        geom_boxplot() + ylim(0,0.45) + ylab("Prop Male") +
        theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


      grace250_aov <- aov(propmale ~ Clone, data = grace250spring)


### Total total graph?

totalCumIndmomwEpp <- week7ACtotalind + week7ACmomPEnormWeek + week7ACMomwEppnormWeek + week7ACPropMalenormWeek +
      week7ACEmbWeek + week7ACmomPEnormCumInd + week7ACMomwEppnormCumInd + week7ACPropMalenormCumInd +
      Methyl + Meso + OneLiter + Grace250 +
      plot_layout(nrow = 3, byrow = TRUE, guides = 'collect') + plot_annotation(tag_levels = 'A')

  totalCumIndmomwEppMesoOnly <- week7ACtotalind + week7ACmomPEnormWeek + week7ACMomwEppnormWeek + week7ACPropMalenormWeek +
      week7ACEmbWeek + week7ACmomPEnormCumInd + week7ACMomwEppnormCumInd + week7ACPropMalenormCumInd +
      plot_layout(nrow = 2, byrow = TRUE, guides = 'collect') + plot_annotation(tag_levels = 'A')

  totalCumIndmomwEppOther <-  Meso + OneLiter + Grace250 + Methyl +
      plot_layout(nrow = 1, byrow = TRUE, guides = 'collect') +
      plot_annotation(tag_levels = list(c("I", "J", "K", "L")))

  totaltotal <- totalCumIndmomwEppMesoOnly / totalCumIndmomwEppOther + plot_layout(heights = c(2, 1)) +
    plot_annotation(tag_levels = 'A')

  newtotal <- week7ACtotalind + week7ACmomPEnormWeek + week7ACMomwEppnormWeek + week7ACPropMalenormWeek +
    week7ACEmbWeek + Methyl + plot_layout(nrow = 2, byrow = TRUE, guides = 'collect') +
    plot_annotation(tag_levels = 'A')

  newtotalB <- week7ACtotalind + week7ACmomPEnormWeek + week7ACMomwEppnormWeek + week7ACPropMalenormWeek +
    week7ACEmbWeekB + Methyl + plot_layout(nrow = 2, byrow = TRUE, guides = 'collect') +
    plot_annotation(tag_levels = 'A')


  ggsave(totaltotal, file="Figure3.pdf")

### Other possible figures/panels

  week7ACTotMale <- ggplot(data=meso7AC, aes(x=Week, y=totmale, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


  week7ACEppTot <- ggplot(data=meso7AC, aes(x=Week, y=LooseEppTotal, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACEppTotbyCumPop <- ggplot(data=meso7AC, aes(x=CumInd, y=LooseEppTotal, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACTotMalebyCumPop <- ggplot(data=meso7AC, aes(x=CumInd, y=totmale, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACEmbDividebyTot <- ggplot(data=meso7AC, aes(x=Week, y=EmbDividebyTot, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACEppDividebyTot<- ggplot(data=meso7AC, aes(x=Week, y=EppDividebyTot, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACEppDividebyTotbyCumPop<- ggplot(data=meso7AC, aes(x=CumInd, y=EppDividebyTot, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  week7ACMaleDividebyTotbyCumPop<- ggplot(data=meso7AC, aes(x=CumInd, y=PropTotMale, color=as.factor(SC), group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
