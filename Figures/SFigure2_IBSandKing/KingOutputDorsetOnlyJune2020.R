#!/usr/bin/env Rscript

### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(tidyverse)

### Read in KING output file
  kinship <- fread("Dorsetking.kin")
  D8PO <- fread("./../D8DBunkParentOffspringRelationships.csv")
  D8POAxCF1 <- D8PO[ParentASC=="A" & ParentBSC=="C"]
  D8POAxCF1subA <- data.table(ID1=D8POAxCF1$FocalIndividual, ID1AxCF1=1)
  D8POAxCF1subB <- data.table(ID2=D8POAxCF1$FocalIndividual, ID2AxCF1=1)

  D8POSubA <- data.table(ID1=D8PO$FocalIndividual, ID1ParentA=D8PO$ParentA, ID1ParentASC=D8PO$ParentASC,
    ID1ParentB=D8PO$ParentB, ID1ParentBSC=D8PO$ParentBSC)
  D8POSubB <- data.table(ID2=D8PO$FocalIndividual, ID2ParentA=D8PO$ParentA, ID2ParentASC=D8PO$ParentASC,
    ID2ParentB=D8PO$ParentB, ID2ParentBSC=D8PO$ParentBSC)

### Read in SC file
  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  sc$population <- str_replace(sc$population, "Dcat", "DCat")
  sc <- sc[Species=="pulex" & Nonindependent==0]
  sc <- sc[Species=="pulex" & Nonindependent==0 & LabGenerated==0]
  #sc <- sc[population=="D8" | population=="DCat" | population=="DBunk"]
  #sc <- sc[population=="D8"]
  #sc <- sc[population=="D8" | population=="DBunk"]

### Add SC info to king output
  scsubA <- data.table(ID1=sc$clone, SCA=sc$SC, medrdA=sc$medrd)
  scsubB <- data.table(ID2=sc$clone, SCB=sc$SC, medrdB=sc$medrd)
  setkey(kinship, ID2)
  setkey(scsubB, ID2)
  tmpm <- merge(kinship, scsubB)
  setkey(tmpm, ID1)
  setkey(scsubA, ID1)
  kinshipsc <- merge(tmpm, scsubA)

### Add in year and population info
  temp <- unlist(strsplit(as.character(kinshipsc$ID1), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  kinshipsc$populationA <- matdat$V3
  kinshipsc$yearA <- matdat$V2

  temp <- unlist(strsplit(as.character(kinshipsc$ID2), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  kinshipsc$populationB <- matdat$V3
  kinshipsc$yearB <- matdat$V2

  kinshipsc$pondcompare <- paste(kinshipsc$populationA, kinshipsc$populationB, sep="_")
  kinshipsc$SCcompare <- paste(kinshipsc$SCA, kinshipsc$SCB, sep="_")
  kinshipsc$clone <- ifelse(kinshipsc$SCA!="OO" & kinshipsc$SCB!="OO" &
    kinshipsc$SCA!="G" & kinshipsc$SCB!="G" &
    kinshipsc$SCA!="selfedC" & kinshipsc$SCB!="selfedC" &
    kinshipsc$SCA!="AxCF1" & kinshipsc$SCB!="AxCF1" &
    kinshipsc$SCA==kinshipsc$SCB, 1, 0)

  setkey(kinshipsc, ID2)
  setkey(D8POAxCF1subB, ID2)
  kinshipscm <- merge(kinshipsc, D8POAxCF1subB, all.x=TRUE)

  setkey(kinshipscm, ID1)
  setkey(D8POAxCF1subA, ID1)
  kinshipscm2 <- merge(kinshipscm, D8POAxCF1subA, all.x=TRUE)

  setkey(kinshipsc, ID2)
  setkey(D8POSubB, ID2)
  kinshipscmallPO <- merge(kinshipsc, D8POSubB, all.x=TRUE)

  setkey(kinshipscmallPO, ID1)
  setkey(D8POSubA, ID1)
  kinshipscmallPO2 <- merge(kinshipscmallPO, D8POSubA, all.x=TRUE)


  kinshipscm2[is.na(ID1AxCF1),ID1AxCF1:=0]
  kinshipscm2[is.na(ID2AxCF1),ID2AxCF1:=0]

  kinshipscmallPO2$PO <- ifelse(kinshipscmallPO2$ID1ParentA==kinshipscmallPO2$ID2 | kinshipscmallPO2$ID1ParentB==kinshipscm2$ID2 |
    kinshipscmallPO2$ID1ParentASC==kinshipscmallPO2$SCB & kinshipscmallPO2$SCB!="OO" | kinshipscmallPO2$ID1ParentBSC==kinshipscmallPO2$SCB & kinshipscmallPO2$SCB!="OO" |
  kinshipscmallPO2$ID2ParentA==kinshipscmallPO2$ID1 | kinshipscmallPO2$ID2ParentB==kinshipscmallPO2$ID1 |
    kinshipscmallPO2$ID2ParentASC==kinshipscmallPO2$SCA & kinshipscmallPO2$SCA!="OO" | kinshipscmallPO2$ID2ParentBSC==kinshipscmallPO2$SCA & kinshipscmallPO2$SCA!="OO", 1, 0)

  kinshipscmallPO2[is.na(PO),PO:=0]

  kinshipscmallPO2$ID1AxCF1 <- ifelse(kinshipscmallPO2$ID1ParentASC=="C" & kinshipscmallPO2$ID1ParentBSC=="A" |
    kinshipscmallPO2$ID1ParentASC=="A" & kinshipscmallPO2$ID1ParentBSC=="C", 1, 0)

  kinshipscmallPO2$ID2AxCF1 <- ifelse(kinshipscmallPO2$ID2ParentASC=="C" & kinshipscmallPO2$ID2ParentBSC=="A" |
    kinshipscmallPO2$ID2ParentASC=="A" & kinshipscmallPO2$ID2ParentBSC=="C", 1, 0)

  kinshipscmallPO2[is.na(ID1AxCF1),ID1AxCF1:=0]
  kinshipscmallPO2[is.na(ID2AxCF1),ID2AxCF1:=0]

  kinshipscmallPO2[is.na(ID1ParentBSC),ID1ParentBSC:="unsampled"]
  kinshipscmallPO2[is.na(ID2ParentBSC),ID2ParentBSC:="unsampled"]


  kinshipscmallPO2$siblings <- ifelse(kinshipscmallPO2$ID1ParentASC==kinshipscmallPO2$ID2ParentASC & kinshipscmallPO2$ID1ParentASC!="OO" &
  kinshipscmallPO2$ID2ParentB!="unsampled" & kinshipscmallPO2$ID1ParentBSC==kinshipscmallPO2$ID2ParentBSC & kinshipscmallPO2$ID1ParentBSC!="OO" |
    kinshipscmallPO2$ID1ParentASC==kinshipscmallPO2$ID2ParentBSC & kinshipscmallPO2$ID1ParentASC!="OO" & kinshipscmallPO2$ID2ParentB!="unsampled" &
    kinshipscmallPO2$ID1ParentBSC==kinshipscmallPO2$ID2ParentASC & kinshipscmallPO2$ID1ParentBSC!="OO" |
    kinshipscmallPO2$ID1ParentA==kinshipscmallPO2$ID2ParentA & kinshipscmallPO2$ID2ParentB!="unsampled" &
    kinshipscmallPO2$ID1ParentB==kinshipscmallPO2$ID2ParentB |
    kinshipscmallPO2$ID1ParentA==kinshipscmallPO2$ID2ParentB & kinshipscmallPO2$ID2ParentB!="unsampled" &
    kinshipscmallPO2$ID1ParentB==kinshipscmallPO2$ID2ParentA,
    1, 0)

  kinshipscmallPO2[is.na(siblings),siblings:=0]


  kinshipscmallPO2$type <- ifelse(kinshipscmallPO2$clone==1, "clone", ifelse(kinshipscmallPO2$clone==0 & kinshipscmallPO2$siblings==1,
    "fullsiblings", ifelse(kinshipscmallPO2$PO==1, "ParentOffspring", ifelse(kinshipscmallPO2$SCA=="poC" &
    kinshipscmallPO2$SCB=="C" | kinshipscmallPO2$SCA=="C" & kinshipscmallPO2$SCB=="poC" | kinshipscmallPO2$SCA=="poH" &
    kinshipscmallPO2$SCB=="H" | kinshipscmallPO2$SCA=="H" & kinshipscmallPO2$SCB=="poH" | kinshipscmallPO2$SCA=="poW" &
    kinshipscmallPO2$SCB=="W" | kinshipscmallPO2$SCA=="W" & kinshipscmallPO2$SCB=="poW" | kinshipscmallPO2$SCA=="B" &
    kinshipscmallPO2$SCB=="poB" | kinshipscmallPO2$SCA=="poB" & kinshipscmallPO2$SCB=="B", "selfedPO", ifelse(kinshipscmallPO2$SCA=="A" &
    kinshipscmallPO2$SCB=="C" | kinshipscmallPO2$SCA=="C" & kinshipscmallPO2$SCB=="A", "AvsC", "other")
    ))))

  kinshipscmallPO2$type <- factor(kinshipscmallPO2$type, levels=c("AvsC", "fullsiblings",
    "clone", "other", "ParentOffspring", "selfedPO"))

  kinshipscmallPO2$toplot <- ifelse(kinshipscmallPO2$type=="other", "other", "plot")

  kinshipscmallPO2$pondcompareB <- ifelse(kinshipscmallPO2$pondcompare=="DBunk_D8", "D8_DBunk", ifelse(
    kinshipscmallPO2$pondcompare=="DCat_D8", "D8_DCat", ifelse(kinshipscmallPO2$pondcompare=="DCat_DBunk",
    "DBunk_DCat", kinshipscmallPO2$pondcompare)
  ))

  kinshipscmallPO2$pondcompareB <- factor(kinshipscmallPO2$pondcompareB, levels=c("D8_D8", "DBunk_DBunk", "DCat_DCat",
    "D8_DBunk", "D8_DCat", "DBunk_DCat"))

  ggplot(data=kinshipscmallPO2[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=type)) + geom_point() +
    geom_point(data = subset(kinshipscmallPO2[medrdA > 9 & medrdB > 9], toplot == "plot")) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  ggplot(data=kinshipscmallPO2[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=type)) + geom_point() +
    geom_point(data = subset(kinshipscmallPO2[medrdA > 9 & medrdB > 9], toplot == "plot")) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(aes(x=0.1336919, y=-0.1272566, color="AvsC")) +
    facet_wrap(~pondcompareB)


    ggplot(data=kinshipscmallPO2[Kinship > 0.2 & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(type))) + geom_point() +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  #kinshipscm2$regionA <- ifelse(kinshipscm2$populationA=="DCat" | kinshipscm2$populationA=="D8" |
  #  kinshipscm2$populationA=="DBunk", "Kilwood", ifelse(kinshipscm2$populationA=="D10", "D10", "W"))
  #kinshipscm2$regionB <- ifelse(kinshipscm2$populationB=="DCat" | kinshipscm2$populationB=="D8" |
  #  kinshipscm2$populationB=="DBunk", "Kilwood", ifelse(kinshipscm2$populationB=="D10", "D10", "W"))


  kinshipscm2$type <- ifelse(kinshipscm2$clone==1, "clone", ifelse(kinshipscm2$clone==0 & kinshipscm2$ID2AxCF1==1 &
    kinshipscm2$ID1AxCF1==1, "AxCF1fullsibs", ifelse(kinshipscm2$ID2AxCF1==1 & kinshipscm2$SCA=="A" |
    kinshipscm2$ID2AxCF1==1 & kinshipscm2$SCA=="C" | kinshipscm2$ID1AxCF1==1 & kinshipscm2$SCB=="A" |
    kinshipscm2$ID1AxCF1==1 & kinshipscm2$SCB=="C", "AxCF1vsACPO", ifelse(kinshipscm2$SCA=="poC" &
    kinshipscm2$SCB=="C" | kinshipscm2$SCA=="C" & kinshipscm2$SCB=="poC" | kinshipscm2$SCA=="poH" &
    kinshipscm2$SCB=="H" | kinshipscm2$SCA=="H" & kinshipscm2$SCB=="poH" | kinshipscm2$SCA=="poW" &
    kinshipscm2$SCB=="W" | kinshipscm2$SCA=="W" & kinshipscm2$SCB=="poW", "selfedPO", ifelse(
    kinshipscm2$SCA=="A" & kinshipscm2$SCB=="C" | kinshipscm2$SCA=="C" & kinshipscm2$SCB=="A", "AvsC",
    ifelse(kinshipscm2$regionA=="Kilwood" & kinshipscm2$regionB=="D10" | kinshipscm2$regionA=="D10" &
    kinshipscm2$regionB=="Kilwood", "KilwoodvsD10", ifelse(kinshipscm2$regionA=="Kilwood" &
    kinshipscm2$regionB=="W" | kinshipscm2$regionA=="W" & kinshipscm2$regionB=="Kilwood", "KilwoodvsW",
    "other"
    )))))))

  kinshipscm2$type <- ifelse(kinshipscm2$clone==1, "clone", ifelse(kinshipscm2$clone==0 & kinshipscm2$ID2AxCF1==1 &
    kinshipscm2$ID1AxCF1==1, "AxCF1fullsibs", ifelse(kinshipscm2$ID2AxCF1==1 & kinshipscm2$SCA=="A" |
    kinshipscm2$ID2AxCF1==1 & kinshipscm2$SCA=="C" | kinshipscm2$ID1AxCF1==1 & kinshipscm2$SCB=="A" |
    kinshipscm2$ID1AxCF1==1 & kinshipscm2$SCB=="C", "AxCF1vsACPO", ifelse(kinshipscm2$SCA=="poC" &
    kinshipscm2$SCB=="C" | kinshipscm2$SCA=="C" & kinshipscm2$SCB=="poC" | kinshipscm2$SCA=="poH" &
    kinshipscm2$SCB=="H" | kinshipscm2$SCA=="H" & kinshipscm2$SCB=="poH" | kinshipscm2$SCA=="poW" &
    kinshipscm2$SCB=="W" | kinshipscm2$SCA=="W" & kinshipscm2$SCB=="poW", "selfedPO", ifelse(
    kinshipscm2$SCA=="A" & kinshipscm2$SCB=="C" | kinshipscm2$SCA=="C" & kinshipscm2$SCB=="A", "AvsC", "other"
    )))))



  kinshipscm2$type <- factor(kinshipscm2$type, levels=c("clone", "selfedPO", "AxCF1vsACPO", "AxCF1fullsibs",
    "other", "AvsC", "KilwoodvsD10", "KilwoodvsW"))

  kinshipscm2$type <- factor(kinshipscm2$type, levels=c("clone", "selfedPO", "AxCF1vsACPO", "AxCF1fullsibs",
    "other", "AvsC"))


  kinshipscm2Kilwoodonly <- kinshipscm2[populationA!="D10" & populationB!="D10" & populationA!="W1" &
    populationB!="W1" & populationA!="W6" & populationB!="W6"]

  ggplot(data=kinshipscm2[Kinship > 0 & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=type)) + geom_point() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  ggplot(data=kinshipscm2[pondcompare=="D8_D8" & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship,
    color=factor(type, labels=c("Clonal", "Selfed parent-offspring",
    "Outcrossed parent-offspring", "Full-siblings", "unknown", "A vs C")))) + geom_point() + labs(color="Type") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  ggplot(data=kinshipscm2[medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship,
    color=factor(type, labels=c("Clonal", "Selfed parent-offspring",
    "Outcrossed parent-offspring", "Full-siblings", "unknown", "A vs C")))) + geom_point() + labs(color="Type") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  ggplot(data=kinshipscm2[type=="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship)) + geom_point() +
    geom_point(data=kinshipscm2[type!="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship,
    color=factor(type, labels=c("Clonal", "Selfed parent-offspring",
    "Outcrossed parent-offspring", "Full-siblings", "A vs C", "Kilwood vs D10", "Kilwood vs W")))) + labs(color="Type") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  ggplot(data=kinshipscm2Kilwoodonly[type=="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship)) + geom_point() +
    geom_point(data=kinshipscm2Kilwoodonly[type!="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship,
    color=factor(type, labels=c("Clonal", "Selfed parent-offspring",
    "Outcrossed parent-offspring", "Full-siblings", "A vs C")))) + labs(color="Type") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  ggplot(data=kinshipscm2[type=="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship)) + geom_point() +
    geom_point(data=kinshipscm2[type!="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship,
    color=factor(type, labels=c("Clonal", "Selfed parent-offspring",
    "Outcrossed parent-offspring", "Full-siblings", "A vs C")))) + labs(color="Type") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


  mean(kinshipscm2$Kinship[kinshipscm2$type=="AvsC"])
  #0.1272566
  kinshipscm2[Kinship<0.1272566 & type!="AvsC"]
  #72910
  kinshipscm2[Kinship>0.1272566 & type!="AvsC"]
  #23644
  72910/(72910+23644)
  #0.7551215

  mean(kinshipscm2$IBS0[kinshipscm2$type=="AvsC"])
  #0.1336919
  kinshipscm2[IBS0<0.1336919 & type!="AvsC"]
  #91522
  kinshipscm2[IBS0>0.1336919 & type!="AvsC"]
  #5032
  91522/(91522+5032)
  #0.9478841

  kinshipscm2[]



  ggplot(data=kinshipscm2[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship,
    color=factor(type))) + geom_point() + labs(color="Type") +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


    test + scale_fill_discrete(
    breaks = c("clone", "selfedPO", "AxCF1vsACPO", "AxCF1fullsibs", "AvsC", "other"),
    labels = c("Clonal", "Selfed parent-offspring", "Outcrossed parent-offpsring", "Full-sibling", "A vs C", "Unassigned")
    )


  kinshipsc$type <- ifelse(kinshipsc$clone==1, "clone", ifelse())


  D8tmp <- kinshipsc[populationB!="DBunk" & populationB!="DBunk" & populationA!="DCat" & populationB!="DCat"]
  D8tmp$type <- ifelse(D8tmp$SCA=="C" & D8tmp$SCB=="poC" | D8tmp$SCA=="poC" & D8tmp$SCB=="C" |
    D8tmp$SCA=="poH" & D8tmp$SCB=="H" | D8tmp$SCA=="H" & D8tmp$SCB=="poH" | D8tmp$SCA=="W" & D8tmp$SCB=="poW" |
    D8tmp$SCA=="poW" & D8tmp$SCB=="W", "selfvsparentfield",
    ifelse(D8tmp$ID1=="May_2017_D8_515" & D8tmp$SCB=="selfedC" | D8tmp$ID2=="May_2017_D8_515" & D8tmp$SCA=="selfedC" |
    D8tmp$ID1=="April_2017_D8_222" & D8tmp$SCB=="selfedC" | D8tmp$ID2=="April_2017_D8_222" & D8tmp$SCA=="selfedC",
    "focalCvsselfedC", ifelse(D8tmp$SCA=="C" & D8tmp$SCB=="selfedC" | D8tmp$SCA=="selfedC" & D8tmp$SCB=="C",
    "nonfocalCvsselfedC", ifelse(D8tmp$clone==1, "clone", "other"))))

  ggplot(data=D8tmp[Kinship > 0 & IBS0 < 0.04], aes(x=IBS0, y=Kinship, color=type)) + geom_point()

  D8tmp$avgmedrd <- (D8tmp$medrdA + D8tmp$medrdB)/2

  ggplot(data=D8tmp[type=="focalCvsselfedC" & N_SNP > 10000 | type=="nonfocalCvsselfedC" & N_SNP > 10000],
    aes(x=N_SNP, y=Kinship, color=type)) + geom_point()

    ggplot(data=D8tmp[type=="focalCvsselfedC" & N_SNP > 10000 | type=="nonfocalCvsselfedC" & N_SNP > 10000],
      aes(x=avgmedrd, y=Kinship, color=type)) + geom_point()

      ggplot(data=D8tmp[type=="focalCvsselfedC" & N_SNP > 10000 | type=="nonfocalCvsselfedC" & N_SNP > 10000],
        aes(x=N_SNP, y=IBS0, color=type)) + geom_point()


  CvsselfedC <- kinshipsc[SCA=="C" & SCB=="selfedC" | SCA=="selfedC" & SCB=="C" |
  CvsselfedC$focalC <- ifelse(CvsselfedC$ID1=="May_2017_D8_515" | CvsselfedC$ID2=="May_2017_D8_515" |
    CvsselfedC$ID1=="April_2017_D8_222" | CvsselfedC$ID2=="April_2017_D8_222", 1, 0)

    ggplot(data=CvsselfedC, aes(x=IBS0, y=Kinship, color=as.factor(focalC))) + geom_point()

  ggplot(data=kinshipsc, aes(x=IBD0, y=IBD1, color=as.factor(clone))) + geom_point()
  ggplot(data=kinshipsc[Kinship > 0.2], aes(x=IBD0, y=IBD1, color=as.factor(clone))) + geom_point()

  kinshipsc$pondcompare <- paste(kinshipsc$populationA,kinshipsc$populationB, sep="_")
  kinshipsc$pondcompareB <- ifelse(kinshipsc$pondcompare=="DBunk_D8", "D8_DBunk", ifelse(
    kinshipsc$pondcompare=="DCat_D8", "D8_DCat", ifelse(kinshipsc$pondcompare=="DCat_DBunk",
    "DBunk_DCat", kinshipsc$pondcompare)
  ))

  ggplot(data=kinshipsc[Kinship > 0 & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(clone))) +
    geom_point() + facet_wrap(~pondcompareB)

    ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(clone))) +
      geom_point() + facet_wrap(~pondcompareB)

  load("./../mclonecountswide_ACfixedSNPs_20200417.Rdata")
  mclonessubID1 <- data.table(ID1=mclonecountswide$clone, ID1_ACF1=mclonecountswide$ACF1hybrid)
  mclonessubID2 <- data.table(ID2=mclonecountswide$clone, ID2_ACF1=mclonecountswide$ACF1hybrid)
  setkey(mclonessubID2, ID2)
  setkey(kinshipsc, ID2)
  mkinshipsc <- merge(kinshipsc, mclonessubID2)
  setkey(mclonessubID1, ID1)
  setkey(mkinshipsc, ID1)
  m2kinshipsc <- merge(mkinshipsc, mclonessubID1)
  m2kinshipsc$ACF1comp <- ifelse(m2kinshipsc$ID2_ACF1=="1" & m2kinshipsc$ID1_ACF1=="1" & m2kinshipsc$clone==0, 1, 0)
  m2kinshipsc$F1vsACparent <- ifelse(m2kinshipsc$SCA=="A" & m2kinshipsc$ID2_ACF1==1 |
    m2kinshipsc$SCA=="C" & m2kinshipsc$ID2_ACF1==1 | m2kinshipsc$SCB=="A" & m2kinshipsc$ID1_ACF1==1 |
    m2kinshipsc$SCB=="C" & m2kinshipsc$ID1_ACF1==1, 1, 0)

  ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9], aes(x=IBD1, y=IBD2, color=as.factor(ACF1comp))) +
    geom_point() + facet_wrap(~pondcompareB)
  ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(F1vsACparent))) +
    geom_point() + facet_wrap(~pondcompareB)

  m2kinshipsc$ACrelation <- ifelse(m2kinshipsc$ACF1comp==1, "ACF1vsACF1", ifelse(
    m2kinshipsc$F1vsACparent==1, "AorCvsACF1", ifelse(m2kinshipsc$SCcompare=="A_C" |
    m2kinshipsc$SCcompare=="C_A", "AvsC", 0)))

    m2kinshipsc$ACrelation <- ifelse(m2kinshipsc$ACF1comp==1, "ACF1vsACF1", ifelse(
      m2kinshipsc$F1vsACparent==1, "AorCvsACF1", ifelse(m2kinshipsc$SCcompare=="A_C" |
      m2kinshipsc$SCcompare=="C_A", "AvsC", ifelse(m2kinshipsc$SCcompare=="C_C" |
      m2kinshipsc$SCcompare=="A_A", "AvsAorCvsC", "neither"))))


  ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9], aes(x=IBD1, y=IBD2, color=as.factor(ACrelation))) +
    geom_point() + facet_wrap(~pondcompareB)

  ggplot(data=m2kinshipsc[populationA=="D8" & populationB=="D8" & medrdA>9 & medrdB>9], aes(x=IBS0, y=Kinship, color=as.factor(ACrelation))) +
    geom_point()

  ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
    aes(x=IBS0, y=Kinship, color=as.factor(ACrelation))) +
    geom_point()

    ggplot(data=m2kinshipsc[IBS0 < 0.025 & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
      aes(x=IBS0, y=Kinship, color=as.factor(ACrelation))) +
      geom_point()

      ggplot(data=m2kinshipsc[IBS0 < 0.025 & medrdA > 5 & medrdB > 5 & populationA=="D8" & populationB=="D8"],
        aes(x=IBS0, y=Kinship, color=as.factor(ACrelation))) +
        geom_point()

  m2kinshipsc$ACrelationB <- ifelse(m2kinshipsc$SCcompare=="H_C" |
    m2kinshipsc$SCcompare=="C_H", "CvsH", ifelse(m2kinshipsc$SCcompare=="W_C" |
    m2kinshipsc$SCcompare=="C_W", "CvsW", ifelse(m2kinshipsc$SCcompare=="W_H" |
    m2kinshipsc$SCcompare=="H_W", "HvsW", m2kinshipsc$ACrelation)))

    ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
      aes(x=IBD1, y=IBD2, color=as.factor(ACrelationB))) +
      geom_point() + facet_wrap(~ACrelationB)

      ggplot(data=m2kinshipsc[Kinship > 0.2 & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
        aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
        geom_point()


        ggplot(data=m2kinshipsc[IBS0 < 0.02 & Kinship > 0.2 & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
          aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
          geom_point()

        ggplot(data=m2kinshipsc[IBS0 < 0.03 & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
          aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
          geom_point()

  m2kinshipsc$yearcompare <- paste(m2kinshipsc$yearA, m2kinshipsc$yearB, sep="_")
  m2kinshipsc$yearcompareB <- ifelse(m2kinshipsc$yearcompare=="2017_2016", "2016_2017", ifelse(
    m2kinshipsc$yearcompare=="2018_2016", "2016_2018", ifelse(
    m2kinshipsc$yearcompare=="2019_2016", "2016_2019", ifelse(
    m2kinshipsc$yearcompare=="2018_2017", "2017_2018", ifelse(
    m2kinshipsc$yearcompare=="2019_2017", "2017_2019", ifelse(
    m2kinshipsc$yearcompare=="2019_2018", "2018_2019", m2kinshipsc$yearcompare
      ))))))

        ggplot(data=m2kinshipsc[yearA==yearB & medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
          aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
          geom_point() + facet_wrap(~yearcompare)

          ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9 & populationA=="D8" & populationB=="D8"],
            aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
            geom_point() + facet_wrap(~yearcompareB)

            ggplot(data=m2kinshipsc[medrdA > 9 & medrdB > 9 & populationA=="DBunk" & populationB=="DBunk"],
              aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
              geom_point() + facet_wrap(~yearcompareB)

              ggplot(data=m2kinshipsc[medrdA > 5 & medrdB > 5 & populationA=="DCat" & populationB=="DCat"],
                aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
                geom_point() + facet_wrap(~yearcompareB)

                ggplot(data=m2kinshipsc[medrdA > 5 & medrdB > 5 & populationA=="D8" & populationB=="D8"],
                  aes(x=IBS0, y=Kinship, color=as.factor(ACrelationB))) +
                  geom_point() + facet_wrap(~yearcompareB)


  PO1 <- m2kinshipsc[Kinship > 0.36 & Kinship < 0.4 & IBS0 < 0.0025]
  PO2 <- m2kinshipsc[Kinship > 0.325 & Kinship < 0.36 & IBS0 < 0.0075]
  PO3 <- m2kinshipsc[Kinship > 0.24 & Kinship < 0.325 & IBS0 < 0.01]
  LikelyClones <- kinshipsc[Kinship > 0.375 & IBS0 < 0.00025]
  LikelyClonesB <- kinshipsc[Kinship > 0.35 & IBS0 < 0.00025]
  LikelyClonesMissed <- LikelyClonesB[clone!=1 & medrdA > 4 & medrdB > 4]

  PO <- rbind(PO1, PO2, PO3)


  PO1 <- kinshipsc[Kinship > 0.3 & Kinship < 0.375 & IBS0 < 0.0019]
  PO2 <- kinshipsc[Kinship > 0.25 & Kinship < 0.3 & IBS0 < 0.0020]
  PO3 <- kinshipsc[Kinship > 0.2 & Kinship < 0.25 & IBS0 < 0.0032]

  kinshipscwinpond <- kinshipsc[populationA==populationB]

  ggplot(data=kinshipscwinpond[Kinship > 0 & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(clone))) +
    geom_point() + facet_wrap(~populationA)

  kinshipscwinpond$AvsC <- ifelse(kinshipscwinpond$SCA=="A" & kinshipscwinpond$SCB=="C" |
    kinshipscwinpond$SCA=="C" & kinshipscwinpond$SCB=="A", 1, ifelse(kinshipscwinpond$ID2=="Spring_2016_D8_8.23"
    & kinshipscwinpond$ID1=="April_2017_D8_214", 2, 0))

  ggplot(data=kinshipscwinpond[Kinship > 0 & medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(AvsC))) +
    geom_point() + facet_wrap(~populationA)

  ggplot(data=kinshipscwinpond[medrdA > 9 & medrdB > 9], aes(x=IBS0, y=Kinship, color=as.factor(AvsC))) +
    geom_point() + facet_wrap(~populationA)

    ggplot(data=kinshipscwinpond[medrdA > 9 & medrdB > 9 & populationA=="D8"], aes(x=IBS0, y=Kinship, color=as.factor(AvsC))) +
      geom_point() + facet_wrap(~AvsC)

  kinshipsc$betweenF1 <- ifelse(kinshipsc$SCcompare=="O_Q" | kinshipsc$SCcompare=="Q_O", 1, ifelse(
    kinshipsc$SCcompare=="O_S" | kinshipsc$SCcompare=="S_O", 2, ifelse(
    kinshipsc$SCcompare=="Q_S" | kinshipsc$SCcompare=="S_Q", 3, ifelse(
    kinshipsc$SCcompare=="O_T" | kinshipsc$SCcompare=="T_O", 4, ifelse(
    kinshipsc$SCcompare=="S_T" | kinshipsc$SCcompare=="T_S", 5, ifelse(
    kinshipsc$SCcompare=="Q_T" | kinshipsc$SCcompare=="T_Q", 6, ifelse(
    kinshipsc$SCcompare=="O_U" | kinshipsc$SCcompare=="U_O", 7, ifelse(
    kinshipsc$SCcompare=="S_U" | kinshipsc$SCcompare=="U_S", 8, ifelse(
    kinshipsc$SCcompare=="Q_U" | kinshipsc$SCcompare=="U_Q", 9, ifelse(
    kinshipsc$SCcompare=="T_U" | kinshipsc$SCcompare=="U_T", 10, ifelse(
    kinshipsc$SCcompare=="O_X" | kinshipsc$SCcompare=="X_O", 11, ifelse(
    kinshipsc$SCcompare=="S_X" | kinshipsc$SCcompare=="X_S", 12, ifelse(
    kinshipsc$SCcompare=="Q_X" | kinshipsc$SCcompare=="X_Q", 13, ifelse(
    kinshipsc$SCcompare=="T_X" | kinshipsc$SCcompare=="X_T", 14, ifelse(
    kinshipsc$SCcompare=="U_X" | kinshipsc$SCcompare=="X_U", 15, 0
  )))))))))))))))

  ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9 & Kinship>0], aes(x=IBS0, y=Kinship, color=as.factor(betweenF1))) +
    geom_point() + facet_wrap(~betweenF1)

    ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9 & Kinship>0], aes(x=IBS0, y=Kinship, color=as.factor(betweenF1))) +
      geom_point()

    kinshipsc$betweenF1B <- ifelse(kinshipsc$betweenF1 > 0, 1, 0)
    ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9 & Kinship>0], aes(x=IBS0, y=Kinship, color=as.factor(betweenF1B))) +
      geom_point()

    ggplot(data=kinshipsc[medrdA > 9 & medrdB > 9 & Kinship>0 & populationA=="D8" & populationB=="D8"],
     aes(x=IBS0, y=Kinship, color=as.factor(betweenF1B))) +
      geom_point()


  kinshipsc <- kinshipsc[populationA!="DOil" & populationA!="Dramp" & populationB!="DOil" & populationB!="Dramp"]
  kinshipsc$withinpond <- ifelse(kinshipsc$populationA==kinshipsc$populationB, 1, 0)
  ggplot(data=kinshipsc, aes(x=as.factor(pondcompare), y=Kinship, color=as.factor(withinpond))) + geom_beeswarm()

### From looking at a graph of Kinship by IBS0, we can see clear groupings falling out. One of these is the line
### along the bottom that appears to respond to clones (individuals from the same clonal lineage.
  # Spring_2016_D8_8.24 appears to belong to superclone H
  # Spring_2016_D8_8.23 appears to belong to superclone C (original B)
  # April_2017_DCat_10 appears to belong to superclone B

  kinshipsc$clone <- ifelse(kinshipsc$SCA=="C" & kinshipsc$ID1=="Spring_2016_D8_8.23" |
    kinshipsc$SCB=="C" & kinshipsc$ID1=="Spring_2016_D8_8.23" |
    kinshipsc$SCA=="C" & kinshipsc$ID2=="Spring_2016_D8_8.23" |
    kinshipsc$SCB=="C" & kinshipsc$ID2=="Spring_2016_D8_8.23", 2, kinshipsc$clone)

    kinshipsc$clone <- ifelse(kinshipsc$SCA=="B" & kinshipsc$ID1=="April_2017_DCat_10" |
      kinshipsc$SCB=="B" & kinshipsc$ID1=="April_2017_DCat_10" |
      kinshipsc$SCA=="B" & kinshipsc$ID2=="April_2017_DCat_10" |
      kinshipsc$SCB=="B" & kinshipsc$ID2=="April_2017_DCat_10", 3, kinshipsc$clone)

    kinshipsc$clone <- ifelse(kinshipsc$SCA=="H" & kinshipsc$ID1=="Spring_2016_D8_8.24" |
      kinshipsc$SCB=="H" & kinshipsc$ID1=="Spring_2016_D8_8.24" |
      kinshipsc$SCA=="H" & kinshipsc$ID2=="Spring_2016_D8_8.24" |
      kinshipsc$SCB=="H" & kinshipsc$ID2=="Spring_2016_D8_8.24", 4, kinshipsc$clone)



  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.19], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()
  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.3], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()

  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.3], aes(x=Kinship, y=IBS0, color=as.factor(clone))) +
    geom_point() + facet_wrap(~clone)

  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.3], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()
  ggplot(data=kinshipsc[N_SNP > 451535 & Kinship > 0.19], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()
  ggplot(data=kinshipsc[N_SNP > 490000 & Kinship > 0.19], aes(x=Kinship, y=IBS0, color=as.factor(clone))) + geom_point()



  msc$pondcompare <- str_replace(msc$pondcompare, "D8_D10", "D10_D8")
  msc$pondcompare <- str_replace(msc$pondcompare, "DBunk_D10", "D10_DBunk")
  msc$pondcompare <- str_replace(msc$pondcompare, "DCat_D10", "D10_DCat")
  msc$pondcompare <- str_replace(msc$pondcompare, "DBunk_D8", "D8_DBunk")
  msc$pondcompare <- str_replace(msc$pondcompare, "DCat_D8", "D8_DCat")
  msc$pondcompare <- str_replace(msc$pondcompare, "W1_D8", "D8_W1")
  msc$pondcompare <- str_replace(msc$pondcompare, "W6_D8", "D8_W6")
  msc$pondcompare <- str_replace(msc$pondcompare, "W1_DBunk", "DBunk_W1")
  msc$pondcompare <- str_replace(msc$pondcompare, "W6_DBunk", "DBunk_W6")
  msc.ag <- msc[,list(meanIBS=mean(Kinship)),
    list(pondcompare, window, chr, start, stop)]
