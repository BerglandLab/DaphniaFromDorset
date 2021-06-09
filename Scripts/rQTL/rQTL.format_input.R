# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(qtl)
  library(data.table)
  #library(ggplot2)
  #library(patchwork)
  library(lme4)

############################
#### Prep the input data ###
############################

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

     male <- mmales[,c("Clone", "Replicate", "Males", "NewTotal", "propmale", "SCB"), with=F]

     ### tack in whether it has been genotyped yet
       #f1s.use <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.onlyPheno.delim")
       #f1s.use <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.allF1s.delim")

       f1s <- sc[AxCF1Hybrid==1][OneLiterPheno==1]$clone
       f1s.use <- data.table(cloneid=f1s)
       setnames(f1s.use, "cloneid", "Clone")
       f1s.use[,geno:=T]

       setkey(mmales, "Clone")
       setkey(f1s.use, "Clone")
       mmales <- merge(mmales, f1s.use, all.x=T, all.y=T)

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
     setkey(f1s.use, "Clone")
     epp <- merge(epp, f1s.use, all.x=T, all.y=T)


### summaries and merging
  ### male
    mmales.ag <- mmales[,list(nClones=length(unique(Clone)),
                              anySeq=any(geno),
                              clone=Clone[geno==T][!is.na(Clone)][1],
                              propmale=sum(Males)/sum(NewTotal),
                              propmalenoneo=sum(Males)/sum(NewTotal-Neos),
                              nMales=mean(Males),
                              nDaps=mean(NewTotal),
                              nNeo=mean(Neos)),
                         list(SCB,
                              population, AxCF1Hybrid)][!is.na(SCB)]


    epp.ag <- epp[,list(nClones=length(unique(Clone)),
                              anySeq=any(geno),
                              fill=mean(fill, na.rm=T),
                              fill.se=sd(fill, na.rm=T)/sqrt(sum(TotalEppB)),
                              epp=mean(TotalEppB, na.rm=T)),
                         list(SCB, population, AxCF1Hybrid)][!is.na(SCB)]



### averages and BLUPs
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


### save
  save(mm, r, mmales, mmales.ag, epp, epp.ag, file="~/F1_pheno.Rdata")
  #load(file="~/F1_pheno.Rdata")


### load data and convert: [M]arkers, [H]eader
  f1s <- fread(file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all_AxC.csv") ### Comes from `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.convert_output.R`
  f1sv <- fread(file= "/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all_AxC.vitPath.csv")
  f1sv[,diplo:=as.numeric(factor(diplo, levels=c("A_m|C_m", "A_p|C_m", "A_m|C_p", "A_p|C_p")))]
  f1s <- f1sv
  f1s[,marker:=paste(chr, id, sep="_")]
  setkey(f1s, marker)

### get one individual per named superclone
  f1s.ind <- data.table(clone=unique(f1s$clone))
  f1s.ind <- merge(f1s.ind, sc, by="clone")
  f1s.ind$SCB <- ifelse(f1s.ind$LabGenerated==0 & f1s.ind$SC!="OO", f1s.ind$SC, f1s.ind$clone)
  f1s.ind.ag <- f1s.ind[OneLiterPheno==1,list(clone=clone[1]), list(SCB)]

  setkey(f1s, clone)
  setkey(f1s.ind.ag, clone)
  f1s.sub <- merge(f1s, f1s.ind.ag)


#### down sample?
  tmp <- data.table(marker=sample(as.character(unique(f1s.sub$marker)), 50000, replace=F))
  setkey(f1s.sub, marker)
  f1s.sub <- f1s.sub[J(tmp)]
  setkey(f1s.sub, id)

  dim(f1s.sub)

#### no?
#  f1s.sub <- f1s

  markers<- dcast(f1s.sub , clone ~ id, value.var=list("diplo"))
  markers[1:5,1:4]


  chrs <- f1s.sub[,list(chr=unique(chr), pos=unique(pos)), id]

  header <- as.data.table(rbind(c("", chrs$chr), c("", chrs$pos)))
  setnames(header, names(header), names(markers))

  markers[1:5,1:4]
  header[1:2,1:4]

  mh <- rbind(header, markers)
  mh[1:5,1:4]

### merge to make rQTL file

  setkey(mm, "clone")
  setkey(mh, "clone")

  mhp <- merge(mm[anySeq==T,c("clone", "SCB",
                              "propmalenoneo", "fill", "epp",
                              "propmalenoneo.ranef", "fill.ranef", "epp.ranef"), with=F],
                mh, all.y=T)

  #setcolorder(mhp, c("clone", "propmalenoneo", "fill", "epp", names(markers)[-1]))

  mhp <- as.matrix(mhp)

  mhp[1:5,1:9]
  mhp[1:2,1:8] <- ""
  mhp[1:5,1:9]

  mhp[1:5,1:9]

  write.csv(mhp, file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/justPheno.rQTL.csv",
            quote=F, row.names=F)






#### test
  f1sv.ag <- f1sv[,list(.N), list(chr, pos, id)]

  ii <- seq(from=1, to=dim(f1sv.ag)[1], by=10)


  out <- foreach(i=ii, .errorhandling="remove")%dopar%{
    print(paste(i, dim(f1sv.ag)[1], sep=" / "))
    tmp <- merge(mm, f1sv[id==f1sv.ag$id[i]], by="clone")

    t1.fill <- anova(lm(fill.ranef~as.factor(geno), tmp))
    t1.male <- anova(lm(propmalenoneo.ranef~as.factor(geno), tmp))


    data.table(id=f1sv.ag$id[i],
               F=c(t1.fill[1,4], t1.male[1,4]),
               p=c(t1.fill[1,5], t1.male[1,5]),
               term=c("fill", "male"))
  }
  out <- rbindlist(out)
  out[,pa:=p.adjust(p, "fdr")]

  save(out, file="~/out.Rdata")



  tmp <- merge(mmales, f1sv[id==981731], by="clone")

  t1.male <- lm(propmalenoneo~as.factor(geno.y), tmp)
  t1.male <- lm(propmalenoneo~as.factor(diplo), tmp)


  summary(t1.male); anova(t1.male)
  save(out, file="~/out.Rdata")

  t1.male <- glmer(propmalenoneo~as.factor(geno.y) + , tmp)

  t1.male <- glm(I(Males/NewTotal)~as.factor(geno.y), family=binomial(), weights=tmp$NewTotal, data=tmp)
  t1.male <- glmer(I(Males/NewTotal)~as.factor(geno.y) + (1|SCB), family=binomial(), weights=tmp$NewTotal, data=tmp)
  t1.male <- lmp(propmalenoneo~as.factor(geno.y), tmp, perm="Exact")
  t1.male <- lmp(propmalenoneo~as.factor(diplo), tmp, perm="Exact")


  library(data.table)
  library(ggplot2)
  load("~/out.Rdata")
  ggplot(data=out, aes(x=id, y=-log10(pa))) + geom_line() + facet_wrap(term~.)
  ggplot(data=out, aes(p)) + geom_histogram() + facet_wrap(term~.)

  out[term=="fill"][which.max(F)]
  tmp <- merge(mm, f1sv[id==out[term=="fill"][which.max(F)]$id], by="clone")

  t1.fill <- (lm(fill.ranef~as.factor(geno), tmp))
  t1.male <- (lm(propmalenoneo.ranef~as.factor(geno), tmp))


### window
  win.bp <- 100000
  step.bp <- 10000
  wins <- foreach(chr.i=unique(f1sv.ag$chr), .combine="rbind")%do%{
    #chr.i=unique(f1sv.ag$chr)[1]
    data.table(chr=chr.i,
              start=seq(from=min(f1sv.ag[chr==chr.i]$pos),
                        to=max(f1sv.ag[chr==chr.i]$pos)-win.bp,
                        by=step.bp),
              stop=seq(from=min(f1sv.ag[chr==chr.i]$pos),
                        to=max(f1sv.ag[chr==chr.i]$pos)-win.bp,
                        by=step.bp)+win.bp)
  }





  f1sv[chr=="Scaffold_9199_HRSCAF_10755"][which.min(abs(pos-6229430))]

  merge(mm, f1sv[pos==6229926], by="clone")

  anova(lm(propmalenoneo.ranef~as.factor(diplo),  merge(mm, f1sv[pos==6229926], by="clone")))
























 #Try glmer
​
   eppB <- epp[fill>=0]
   mmalesB <- mmales[propmale>0]
   #t1 <- glmer(propmale~1+(1|Replicate)+(1|Clone)+(1|SCB), data=mmales, family=(binomial), weights=NewTotal)
   #t2 <- glmer(fill~1+(1|Rep)+(1|Clone), data=eppB, family=binomial(), weights=TotalEppB)
​
   #t1 <- glmer(propmale~1+(1|SCB:Clone:Replicate), data=mmalesB, family=(binomial), weights=NewTotal)
   #t2 <- glmer(propmale~1+(1|Clone:Replicate), data=mmalesB, family=(binomial), weights=NewTotal)
​
   t1 <- glmer(propmale~1+(1|Clone:Replicate) + (1|SCB), data=mmales, family=(binomial), weights=NewTotal)
   t2 <- glmer(propmale~1+(1|Clone:Replicate), data=mmales, family=(binomial), weights=NewTotal)
​
​
   anova(t1, t2)
​
   Data: mmales
   Models:
   t2: propmale ~ 1 + (1 | Clone:Replicate)
   t1: propmale ~ 1 + (1 | Clone:Replicate) + (1 | SCB)
      npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
   t2    2 623.61 629.08 -309.81   619.61
   t1    3 598.62 606.82 -296.31   592.62 26.996  1  2.038e-07 ***
​
​
   malesresid <- data.table(SCB=dimnames(ranef(t1)[[2]])[[1]] ,
                     malesresid=ranef(t1)[[2]][,1] + fixef(t1))
​
   males.ag <- mmales[,list(meanmalesbyclone=mean(propmale)),
     list(Clone, SCB)]
​
   malesbySCB.ag <- males.ag[,list(meanmalesbySCB=mean(meanmalesbyclone)),
     list(SCB)]
​
     popsize.ag <- mmales[,list(meantotal=mean(NewTotal)),
       list(SCB)]
​
​
   setkey(males.ag, SCB)
   setkey(malesresid, SCB)
   m <- merge(males.ag, malesresid)
​
   ggplot(data=m, aes(x=meanmales, y=malesresid)) + geom_point()
​
   setkey(popsize.ag, SCB)
   setkey(malesresid, SCB)
   m2 <- merge(popsize.ag, malesresid)
​
   ggplot(data=m2, aes(x=meantotal, y=malesresid)) + geom_point()
​
​
   setkey(m2pheno, SCB)
   setkey(malesresid, SCB)
   m3 <- merge(m2pheno, malesresid)
​
   ggplot(data=m3, aes(x=embresid, y=malesresid)) + geom_point()
​
​
​
   epp <- fread("EphippiaFinal.csv")
   sc <- fread("/Users/kbkubow/Box Sync/Daphnia/InitialManuscript/SupercloneFiles/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
   sconeliter <- sc[OneLiterPheno==1]
   sconeliter$SCB <- ifelse(sconeliter$LabGenerated==0 & sconeliter$SC!="OO", sconeliter$SC, sconeliter$clone)
   sconeliter$Clone <- sconeliter$clone
​
   setkey(epp, Clone)
   setkey(sconeliter, Clone)
   mepp <- merge(epp, sconeliter)
   epp <- mepp
​
 # Add new variables
   epp$TotalEppB <- epp$TotalEpp*2
   epp$fill <- epp$TotalEmbCorr/epp$TotalEppB
   epp$fullid <- paste(epp$Clone,epp$Rep,sep="_")
   epp$Rep <- as.factor(epp$Rep)
   epp$Clone <- as.factor(epp$Clone)
   epp$Dayssince1stDate <- as.factor(epp$Dayssince1stDate)
   epp$SCB <- as.factor(epp$SCB)
​
   epp.ag <- epp[,list(meaneppbyclone=mean(TotalEpp)),
     list(Clone, SCB)]
​
   eppbySCB.ag <- epp.ag[,list(meaneppbySCB=mean(meaneppbyclone)),
     list(SCB)]
​
   emb.ag <- epp[,list(meanembbyclone=mean(TotalEmbCorr)),
     list(Clone, SCB)]
​
   embbySCB.ag <- emb.ag[,list(meanembbySCB=mean(meanembbyclone)),
     list(SCB)]
​
   setkey(eppbySCB.ag, SCB)
   setkey(embbySCB.ag, SCB)
   mmmeans <- merge(eppbySCB.ag, embbySCB.ag)
   setkey(mmmeans, SCB)
   setkey(malesbySCB.ag, SCB)
   m2means <- merge(mmmeans, malesbySCB.ag, all.x=TRUE)
   save(m2means, file="m2means.Rdata")
​
   t3 <- glmer(fill~1+(1|Dayssince1stDate)+(1|Clone:Rep)+(1|SCB), data=eppB, family=(binomial), weights=TotalEppB)
   t4 <- glmer(fill~1+(1|Dayssince1stDate)+(1|Clone:Rep), data=eppB, family=(binomial), weights=TotalEppB)
​
   anova(t3, t4)
​
   embresid <- data.table(SCB=dimnames(ranef(t3)[[2]])[[1]] ,
                     embresidB=ranef(t3)[[2]][,1] + fixef(t3))
​
   setkey(m3, SCB)
   setkey(embresid, SCB)
   m4 <- merge(m3, embresid)
​
   ggplot(data=m4, aes(x=embresid, embresidB)) + geom_point()
​
​
   t5 <- glmer(TotalEpp~1+(1|Dayssince1stDate)+(1|Clone:Rep)+(1|SCB), data=epp, family=(poisson))
   t6 <- glmer(TotalEpp~1+(1|Dayssince1stDate)+(1|Clone:Rep), data=epp, family=(poisson))
​
   anova(t5, t6)
​
   eppresid <- data.table(SCB=dimnames(ranef(t5)[[2]])[[1]] ,
                     eppresidB=ranef(t5)[[2]][,1] + fixef(t5))
​
   setkey(m4, SCB)
   setkey(eppresid, SCB)
   m5 <- merge(m4, eppresid)
​
   setkey(eppresid, SCB)
   setkey(embresid, SCB)
   setkey(malesresid, SCB)
   mpheno <- merge(eppresid, embresid, all.x=TRUE)
   setkey(mpheno, SCB)
   mphenoupdate <- merge(mpheno, malesresid, all.x=TRUE)
​
   save(mphenoupdate, file="mphenoupdate_20200824.Rdata")
​
   ggplot(data=mphenoupdate, aes(x=malesresid, embresidB)) + geom_point()
​




### set up [P]henotype data
  ### This uses the raw data
    # load(file="~/m3epp.Rdata")
    # phenos <- m3epp.ag[Type=="AxCF1", c("mu.epp", "mu.fill", "Clone"), with=F]

    # phenos[Clone=="D818111", Clone:="April5_2018_D8_18111"]
    # phenos[Clone=="D818106", Clone:="April5_2018_D8_18106"]
    # phenos[Clone=="D818025", Clone:="March20_2018_D8_18025"]
    # phenos[Clone=="D818028", Clone:="March20_2018_D8_18028"]
    # phenos[Clone=="D818030", Clone:="March20_2018_D8_18030"]
    # phenos[Clone=="D818010", Clone:="March20_2018_D8_18010"]

    # setnames(phenos, "Clone", "clone")

  ### This uses the BLUP phenotypes, extracted from the rQTL file that Karen had made
    AxCF1 <- read.cross("csv","","~/AxCF1genoandphenoreconstruct_sub3.csv", genotypes=NULL)
    AxCF1$pheno

    phenos <- as.data.table(AxCF1$pheno)
    setnames(phenos, "CloneID", "clone")

    ### ---> Tack in new phenotype data here

### merge to make rQTL file
  setkey(phenos, "clone")
  setkey(mh, "clone")

  mhp <- merge(phenos, mh, all.y=T)

  setcolorder(mhp, c("eppresid", "embresid", "SCB", "clone", names(markers)[-1]))

  mhp[1:5,1:4]
  dim(mhp)

  write.csv(mhp, file="~/mhp.csv", quote=F, row.names=F, na="")

### Run rQTL
  ### read cross object
    AxCF1 <- read.cross("csv","","~/mhp.csv", crosstype="4way", genotypes=NULL)

  ### marker regression

    mr1 <- scanone(AxCF1, pheno.col=1, method="mr")
    mr2 <- scanone(AxCF1, pheno.col=2, method="mr")

    mr1.dt <- as.data.table(mr1)
    mr2.dt <- as.data.table(mr2)

### Load pooled data
    load("~/peaks.Rdata")
    peaks <- fread("/Users/alanbergland/peaks.csv")

    setnames(peaks, "CHROM", "chr")
    setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))

### plot
    eppresid.plot <- ggplot() +
    #geom_hline(yintercept=summary(perm)[2]) +
    geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
    geom_line(data=mr1.dt, aes(x=pos, y=lod, color=chr), size=1) +
    facet_grid(.~chr, scales="free_x") +
    theme(strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
    ggtitle("epp.resid")


    #setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))
    Gprime.plot <- ggplot(data=gprime, aes(x=pos, y=Gprime, color=chr)) +
    geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
    geom_line(size=.75) +
    facet_grid(.~chr, scales="free_x") +
    theme(strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
    ggtitle("pooledWild")

    embresid.plot <- ggplot() +
    #geom_hline(yintercept=summary(perm)[2]) +
    geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
    geom_line(data=mr2.dt, aes(x=pos, y=lod, color=chr), size=1) +
    facet_grid(.~chr, scales="free_x") +
    theme(strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
    ggtitle("emb.resid")


### final plot
  eppresid.plot / Gprime.plot / embresid.plot


### save
  tar czvf rQTL.inputFiles.tar.gz \
  ~/AxC_F1.csv \
  ~/AxCF1genoandphenoreconstruct_sub3.csv \
  ~/peaks.Rdata \
  /Users/alanbergland/peaks.csv
