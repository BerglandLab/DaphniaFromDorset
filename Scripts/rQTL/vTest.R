library(data.table)

load("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")
male.ag <- male[!is.na(gr),list(propmale=sum(Males)/sum(NewTotal),
                                N=sum(NewTotal),
                                male.mean=mean(Males/NewTotal),
                                male.sd=sd(Males/NewTotal),
                                n=length(Males)),
                  list(clone, gr)]

### parental means and variances
  parA <- glmer(cbind(Males, NewTotal-Males) ~ 1 + (1|Replicate), male[gr=="A"], family=binomial())
  parC <- glmer(cbind(Males, NewTotal-Males) ~ 1 + (1|Replicate), male[gr=="C"], family=binomial())
  AxC <- glmer(cbind(Males, NewTotal-Males) ~ 1 + (1|Replicate) + (1|clone), male[gr=="AxC"], family=binomial())
  #CxC <- glmer(cbind(Males, NewTotal-Males) ~ 1 + (1|Replicate:clone) + (1|clone), male[gr=="CxC"], family=binomial())

  sig.sq.par <- var(c(fixef(parA), fixef(parC)))
  n.parA <- dim(ranef(parA)$Replicate)[1]
  n.parC <- dim(ranef(parC)$Replicate)[1]
  sd.parA <- VarCorr(parA)$Replicate[1,1]
  sd.parC <- VarCorr(parC)$Replicate[1,1]

  sd.AxC <- VarCorr(AxC)$clone[1,1]

  H2 <- (VarCorr(AxC)$clone[1,1] - VarCorr(AxC)$Replicate[1,1]) / VarCorr(AxC)$clone[1,1]
  c <- 2

  v <- (sig.sq.par - sd.parA/n.parA - sd.parC/n.parC) / (sd.AxC * H2 * c)
  v

  n <- (sig.sq.par - sd.parA/n.parA - sd.parC/n.parC) / (8*sd.AxC)
  n

### fill rate
parA <- glmer(cbind(fill*TotalEppB, TotalEppB*(1-fill)) ~ 1 + (1|Rep), epp[gr=="A"], family=binomial())
parC <- glmer(cbind(fill*TotalEppB, TotalEppB*(1-fill)) ~ 1 + (1|Rep), epp[gr=="C"], family=binomial())
AxC <- glmer(cbind(fill*TotalEppB, TotalEppB*(1-fill)) ~ 1 + (1|Rep) + (1|clone), epp[gr=="AxC"], family=binomial())
CxC <- glmer(cbind(fill*TotalEppB, TotalEppB*(1-fill)) ~ 1 + (1|Rep) + (1|clone), epp[gr=="CxC"], family=binomial())

sig.sq.par <- var(c(fixef(parA), fixef(parC)))
n.parA <- dim(ranef(parA)$Rep)[1]
n.parC <- dim(ranef(parC)$Rep)[1]
sd.parA <- VarCorr(parA)$Rep[1,1]^2
sd.parC <- VarCorr(parC)$Rep[1,1]^2

sd.AxC <- VarCorr(AxC)$clone[1,1]^2

H2 <- (VarCorr(AxC)$clone[1,1] - VarCorr(AxC)$Rep[1,1]) / VarCorr(AxC)$clone[1,1]
c <- 2

v <- (sig.sq.par - sd.parA/n.parA - sd.parC/n.parC) / (sd.AxC * H2 * c)
v

n <- (sig.sq.par - se.parA - se.parC) / (8*sd.AxC)


  epp$gr <- ifelse(epp$SC=="selfedA", "A", ifelse(epp$SC=="B", "B", epp$gr))
  epp$counts <- c(1)

  epp.ag <- epp[!is.na(gr),list(fillrate=sum(fill*TotalEppB, na.rm=T)/sum(TotalEppB),
                                N=sum(TotalEppB), meanEMB=mean(TotalEmb),
                                sdEMB=sd(TotalEmb), sampsize=sum(counts),
                                meanEpp=mean(TotalEpp), sdEpp=sd(TotalEpp)),
                  list(clone, gr)]







































#library(qtl)
library(data.table)
library(lme4qtl)
library(doMC)
registerDoMC(10)

############################
#### Prep the input data ###
############################
### `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.convert_output.vitPath.R` does this

########################
### [P]henotype data ###
########################
### Load SC data
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020//Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

## Male rates
  ### load raw male phenotype data
    males <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/MaleCensus.csv")
    #males <- fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rqtl.convertData.R")

  ### some renaming
   sconeliter <- sc[OneLiterPheno==1]
   sconeliter$SCB <- ifelse(sconeliter$LabGenerated==0 & sconeliter$SC!="OO", sconeliter$SC, sconeliter$clone)
   sconeliter$Clone <- sconeliter$clone

  ### tidy
   setkey(males, Clone)
   setkey(sconeliter, Clone)
   mmales <- merge(males, sconeliter)
   mmales$propmale <- mmales$Males/mmales$NewTotal
   mmales$propmalenoneo <- mmales$Males/(mmales$NewTotal-mmales$Neos)

   male <- mmales[,c("Clone", "Replicate", "Males", "NewTotal", "Neos", "propmale", "SC", "SCB", "AxCF1Hybrid", "OneLiterPheno"), with=F]
   setnames(male, "Clone", "clone")


## Epphipia fill rates and production
  ### load raw epphiphial fill data
   epp <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/EphippiaFinal.csv")
   sconeliter <- sc[OneLiterPheno==1]
   sconeliter$SCB <- ifelse(sconeliter$LabGenerated==0 & sconeliter$SC!="OO", sconeliter$SC, sconeliter$clone)
   sconeliter$Clone <- sconeliter$clone

   setkey(epp, Clone)
   setkey(sconeliter, Clone)
   mepp <- merge(epp, sconeliter)
   epp <- mepp

 # Add new variables
   epp$TotalEppB <- epp$TotalEpp*2
   epp$fill <- epp$TotalEmbCorr/epp$TotalEppB
   epp$fullid <- paste(epp$Clone,epp$Rep,sep="_")
   epp$Rep <- as.factor(epp$Rep)
   epp$Clone <- as.factor(epp$Clone)
   epp$Dayssince1stDate <- as.factor(epp$Dayssince1stDate)
   epp$SCB <- as.factor(epp$SCB)

   setkey(epp, "Clone")

### add in cross type
  epp[SCB=="A", gr:="A"]
  epp[SCB=="C", gr:="C"]
  epp[AxCF1Hybrid==1, gr:="AxC"]
  epp[OneLiterPheno==1 & AxCF1Hybrid==0 & SC=="selfedC", gr:="CxC"]

  table(epp$gr)

  male[SCB=="A", gr:="A"]
  male[SCB=="C", gr:="C"]
  male[AxCF1Hybrid==1, gr:="AxC"]
  male[OneLiterPheno==1 & AxCF1Hybrid==0 & SC=="selfedC", gr:="CxC"]

  table(epp$gr)

  table(male$gr)


  epp


























### a few plot
  # scp aob2x@rivanna.hpc.virginia.edu:~/F1_pheno.Rdata ~/.
  # R
  library(ggplot2)
  library(data.table)
  library(patchwork)
  library(cowplot)
  theme_set(theme_cowplot())
  library(lme4)

  load("~/F1_pheno.Rdata")

  ### averages
    epp.ag[SCB=="A", gr:="A"]
    epp.ag[SCB=="C", gr:="C"]
    epp.ag[(grepl("AxB", SCB) | AxCF1Hybrid==1) & anySeq==T, gr:="AxC"]
    epp.ag[grepl("AxB", SCB) & is.na(anySeq), gr:="CxC"]

    setkey(mmales.ag, SCB, population, anySeq)
    setkey(epp.ag, SCB, population, anySeq)


  ### generate BLUPs
  male.mod <- glmer(propmalenoneo ~ 1 + (1|SCB) + (1|Clone:Replicate),
                    data=mmales,
                    family=binomial(),
                    weights=NewTotal-Neos)

  fill.mod <- glmer(fill~1+(1|Dayssince1stDate)+(1|Clone:Rep)+(1|SCB),
              data=epp,
              family=binomial(),
              weights=TotalEppB)

  epp.mod <- glmer(TotalEppB~1+(1|Dayssince1stDate)+(1|Clone:Rep)+(1|SCB),
              data=epp,
              family=poisson())


  r1 <- data.table(SCB=rownames(ranef(male.mod)$SCB),
                  propmalenoneo.ranef =  ranef(male.mod)$SCB[,1])
  r2 <- data.table(SCB=rownames(ranef(fill.mod)$SCB),
                  fill.ranef =  ranef(fill.mod)$SCB[,1])
  r3 <- data.table(SCB=rownames(ranef(epp.mod)$SCB),
             epp.ranef =  ranef(epp.mod)$SCB[,1])

  r <- merge(r1, r2, by="SCB", all.x=T, all.y=T)
  r <- merge(r, r3, by="SCB", all.x=T, all.y=T)


  m <- merge(mmales.ag, epp.ag, all.x=T, all.y=T)

  mm <- merge(m, r, by="SCB")

plot(propmalenoneo~propmalenoneo.ranef, mm)
plot(fill~fill.ranef, mm)
plot(epp~epp.ranef, mm)


  ### weighted means
    fill.rate <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=fill, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05), size=5) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Ephippial Fill Rate")

    num.epp <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=epp, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05), size=5) +
    theme(legend.position = "none")+
    xlab("") +
    ylab("Number of Ephippia")

    propmale <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=propmalenoneo, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05), size=5) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Proportion Male")

  ### ranef
    fill.rate.r <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=fill.ranef, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05), size=5) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Ephippial Fill Rate, BLUP")

    num.epp.r <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=epp.ranef, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05), size=5) +
    theme(legend.position = "none")+
    xlab("") +
    ylab("Number of Ephippia, BLUP")

    propmale.r <- ggplot(data=mm[!is.na(gr)], aes(x=gr, y=propmalenoneo.ranef, color=as.factor(gr))) +
    geom_point(position=position_jitter(width=0.05), size=5) +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Proportion Male, BLUP")


    (fill.rate | num.epp | propmale) /
    (fill.rate.r| num.epp.r | propmale.r)



### princomp
  pr <- prcomp(na.omit(mm[,c("fill", "epp", "propmalenoneo"), with=F]), center=T, scale=T)

  ggplot(data=mm[!is.na(gr)], aes(y=fill, x=propmalenoneo, color=as.factor(gr))) +
  geom_point()
