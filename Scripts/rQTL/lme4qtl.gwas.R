# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

args = commandArgs(trailingOnly=TRUE)
perm <- as.numeric(args[1]) - 1
set <- (args[2]) ### all,  AxC, or CxC

#perm <- 0; set<-"AxC"

### libraries
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

    #save(male, epp, file="/nv/vol186/bergland-lab/alan/phenos_F1.Rdata")

    ### subset down
      if(set=="AxC") {
        epp <- epp[gr=="AxC"]
        male <- male[gr=="AxC"]
      } else if (set=="CxC") {
        epp <- epp[gr=="CxC"]
        male <- male[gr=="CxC"]
      } else if (set=="all") {
        epp <- epp[gr%in%c("AxC", "CxC")]
        male <- male[gr%in%c("AxC", "CxC")]
      }



    #male.full <- glmer(cbind(Males, NewTotal-Males-Neos) ~ 1 + (1|clone) + (1|Replicate),
    #                    family=binomial(),
    #                    data=male[gr=="AxC"])
    #male.red <- glmer(cbind(Males, NewTotal-Males-Neos) ~ 1 + (1|Replicate),
    #                    family=binomial(),
    #                    data=male[gr=="AxC"])
    #anova(male.full, male.red)


#####################
### Genotype data ###
#####################

  ### load genotype data (from DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.convert_output.vitPath.R)
    pio <- fread(file= "/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all.vitPath.pio.csv")
    pio[,phase.fold.geno:=phase.geno]
    pio[phase.fold.geno%in%c(12, 21), phase.fold.geno:=12]

  ### trim to samples used
    setkey(pio, clone)
    pio <- pio[J(unique(c(epp$clone, male$clone)))]

  ### identify unique SNPs (there shoudl be about 5K of them)
    pio.unique <- pio[,list(vec=paste(phase.fold.geno, collapse=".")), list(chr, pos, id)]
    pio.uniq.tab <- table(pio.unique$vec)
    pio.uniq.n <- pio.unique[,list(id=id[1], chr=chr[1], pos=pos[1], .N,
                                  ids=paste(id, collapse=";"),
                                  chrs=paste(chr, collapse=";"),
                                  poss=paste(pos, collapse=";")),
                              list(vec)]

    pio.uniq.n[,set:=set]
    save(pio.uniq.n, file=paste("/scratch/aob2x/daphnia_hwe_sims/lmer4qtl/pio.uniq.set.", set, ".Rdata", sep=""))

    setkey(pio, id, chr, pos)
    setkey(pio.uniq.n, id, chr, pos)
    pio.u <- merge(pio, pio.uniq.n[,-c("vec", "ids"), with=F])


    ((table(pio.u$imputedGeno, pio.u$phase.fold.geno)))
    ((table(pio.u$obs.dosage, pio.u$phase.fold.geno)))


  ### make PIO unique into wide
      pio.u.ag <-   pio.u[,list(phase.fold.geno=phase.fold.geno[1],
                            ng=length(unique(phase.fold.geno))),
                      list(chr, pos, id, clone)]

      pio.u.ag[,d:=as.numeric(factor(phase.fold.geno)) - 2]
      pio.u.ag.ag <- pio.u.ag[,list(maf=mean((d+1)/2)), list(id)]

      #setkey(pio.u.ag, id)
      #pio.u.ag <- pio.u.ag[J(pio.u.ag.ag[maf>.05 & maf<.95]$id)]


      pio.u.ag.w <- dcast(pio.u.ag, clone ~ id, value.var="d")
      pio.u.ag.ag <- pio.u.ag[,list(chr=chr[1], pos=pos[1]), list(id)]


#### using lme4qtl
  lmer.gwas <- foreach(i=pio.u.ag.ag$id, .errorhandling="remove")%dopar%{

    message(which(i==pio.u.ag.ag$id))
    #i<-pio.u.ag.ag$id[9]
    #i=2920958
    #i=984400

    # i= pio.u.ag.ag[chr=="Scaffold_9199_HRSCAF_10755"][which.min(abs(pos-6232363))]$id

    ### make new kinship matrix
      kinship_matrix <- rrBLUP::A.mat(pio.u.ag.w[,-c("clone", i),with=F], n.core = 10)
      dimnames(kinship_matrix)[[1]] <- (pio.u.ag.w$clone)
      dimnames(kinship_matrix)[[2]] <- (pio.u.ag.w$clone)



    geno.tmp <- pio.u[id==i]

    ### male
      #perm <- 0
      if(perm==0) {
        gp.male.tmp <- merge(male, geno.tmp, by="clone")

      } else if(perm>0) {
        set.seed(perm)
        male.ag <- male[,list(.N), list(clone)]
        male.ag[,permClone:=sample(clone)]

        male.perm <- merge(male, male.ag, by="clone")
        setnames(male.perm, c("clone", "permClone"), c("origClone", "clone"))

        gp.male.tmp <- merge(male.perm, geno.tmp, by="clone")

      }

      gp.male.tmp[,d:=as.numeric(as.factor(phase.fold.geno))-2]

      if(set=="AxC" | set=="CxC") {
        male.full <- relmatGlmer(cbind(Males, NewTotal-Males-Neos) ~ d + (1|clone) + (1|Replicate),
                            family=binomial(),
                            data=gp.male.tmp,
                            relmat=list(clone=kinship_matrix))

        male.red <- relmatGlmer(cbind(Males, NewTotal-Males-Neos) ~ 1 + (1|clone) + (1|Replicate),
                            family=binomial(),
                            data=gp.male.tmp,
                            relmat=list(clone=kinship_matrix))
      } else if(set=="all") {
        male.full <- relmatGlmer(cbind(Males, NewTotal-Males-Neos) ~ d + (1|clone) + (1|Replicate),
                            family=binomial(),
                            data=gp.male.tmp,
                            relmat=list(clone=kinship_matrix))

        male.red <- relmatGlmer(cbind(Males, NewTotal-Males-Neos) ~ 1 + (1|clone) + (1|Replicate),
                            family=binomial(),
                            data=gp.male.tmp,
                            relmat=list(clone=kinship_matrix))

      }
      aov.male <- anova(male.full, male.red)



    ### fill rate
      if(perm==0) {
        gp.fill.tmp <- merge(epp, geno.tmp, by="clone")

      } else if (perm>0) {
        set.seed(perm)
        epp.ag <- epp[,list(.N), list(clone)]
        epp.ag[,permClone:=sample(clone)]

        epp.perm <- merge(epp, epp.ag, by="clone")
        setnames(epp.perm, c("clone", "permClone"), c("origClone", "clone"))

        gp.fill.tmp <- merge(epp.perm, geno.tmp, by="clone")

      }

      gp.fill.tmp[,d:=as.numeric(as.factor(phase.fold.geno))-2]

      if(set=="AxC" | set=="CxC") {
        fill.full <- relmatGlmer(cbind(fill*TotalEppB, (1-fill)*TotalEppB) ~ d  + (1|clone) + (1|Rep),
                            family=binomial(),
                            data=gp.fill.tmp,
                            relmat=list(clone=kinship_matrix))

        fill.red <- relmatGlmer(cbind(fill*TotalEppB, (1-fill)*TotalEppB) ~ 1 + (1|clone) + (1|Rep),
                            family=binomial(),
                            data=gp.fill.tmp,
                            relmat=list(clone=kinship_matrix))
      } else if(set=="all") {
        fill.full <- relmatGlmer(cbind(fill*TotalEppB, (1-fill)*TotalEppB) ~ d  + (1|clone) + (1|Rep),
                            family=binomial(),
                            data=gp.fill.tmp,
                            relmat=list(clone=kinship_matrix))

        fill.red <- relmatGlmer(cbind(fill*TotalEppB, (1-fill)*TotalEppB) ~ 1 + (1|clone) + (1|Rep),
                            family=binomial(),
                            data=gp.fill.tmp,
                            relmat=list(clone=kinship_matrix))

      }
      aov.fill <- anova(fill.full, fill.red)

    ### ephpippial number
      if(set=="AxC" | set=="CxC") {
        epp.full <- relmatGlmer(I(TotalEppB/2) ~ d  + (1|clone) + (1|Rep),
                            family=poisson(),
                            data=gp.fill.tmp,
                            relmat=list(clone=kinship_matrix))

        epp.red <- relmatGlmer(I(TotalEppB/2) ~ 1 + (1|clone) + (1|Rep),
                            family=poisson(),
                            data=gp.fill.tmp,
                            relmat=list(clone=kinship_matrix))
      } else if(set=="all") {
        epp.full <- relmatGlmer(I(TotalEppB/2) ~ d  + (1|clone) + (1|Rep),
                            family=poisson(),
                            data=gp.fill.tmp,
                            relmat=list(clone=kinship_matrix))

        epp.red <- relmatGlmer(I(TotalEppB/2) ~ 1 + (1|clone) + (1|Rep),
                            family=poisson(),
                            data=gp.fill.tmp,
                            relmat=list(clone=kinship_matrix))

      }
      aov.epp <- anova(epp.full, epp.red)



    ### output
      averages <- merge(gp.male.tmp[,list(males.mu=sum(Males)/sum(NewTotal-Neos), males.n=sum(NewTotal-Neos)), list(phase.fold.geno, d)],
                        gp.fill.tmp[,list(fill.mu=sum(fill*TotalEppB, na.rm=T)/sum((1-fill)*TotalEppB, na.rm=T), fill.n=sum(TotalEppB, na.rm=T)), list(phase.fold.geno, d)],
                        by="phase.fold.geno")
      averages[,id:=i]

      out <- data.table(id=i,
                  chisq=c(aov.fill[2,6], aov.epp[2,6], aov.male[2,6]),
                  p.aov=c(aov.fill[2,8], aov.epp[2,8], aov.male[2,8]),
                  p.z=c(summary(fill.full)$coef[2,4], summary(epp.full)$coef[2,4], summary(male.full)$coef[2,4]),
                  b.z=c(summary(fill.full)$coef[2,1], summary(epp.full)$coef[2,1], summary(male.full)$coef[2,1]),
                  term=c("fill", "epp", "male"),
                  sing=any(lme4::isSingular(fill.full),
                           lme4::isSingular(fill.red),
                           lme4::isSingular(epp.full),
                           lme4::isSingular(epp.red),
                           lme4::isSingular(male.full),
                           lme4::isSingular(male.red)),
                  set=paste(unique(gp.fill.tmp$gr), collapse="."))
      out[,perm:=perm]
      merge(out, averages, by="id", allow.cartesian=T)

  }
  lmer.gwas <- rbindlist(lmer.gwas)

  lmer.gwas <- merge(lmer.gwas, pio.u.ag.ag, by="id")



  save(lmer.gwas, file=paste("/scratch/aob2x/daphnia_hwe_sims/lmer4qtl/v3.perm", perm, ".set.", set, ".Rdata", sep=""))
