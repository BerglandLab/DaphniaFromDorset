### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(data.table)
  library(SeqArray)
  library(SNPRelate)

############
### data ###
############

### open GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)
  #genofile <- snpgdsOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.gds", allow.duplicate=TRUE)

### which F1s?
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

### load summarized GWAS output
  load("~/lme4qtl_output.AxC.long.SUMMARIZED.Rdata")


########################
### get summary info ###
########################

  ### identify representative clone id per SC for phenotyped strains
    target.sampleId <- epp[,list(.N, sampleid=Clone[which.max(medrd)]), SCB]

  ### get positions of marker SNPs
    target.variantId <-  o.ag.plot[,list(.N, id=min(id), sig=(p.aov<=p.aov.thr)[1],
                                          chr=chr[which.min(id)], pos=pos[which.min(pos)]),
                                    list(chisq, term)]

################
### analysis ###
################

  LDmat <- snpgdsLDMat(genofile,
                      sample.id=target.sampleId$sampleid,
                      snp.id=unique(target.variantId$id),
                      slide=0,
                      method=c("corr"), mat.trim=T,
                      num.thread=1L, with.id=TRUE, verbose=TRUE)

  LDlong <- data.table(cor=expand.grid(LDmat$LD)$Var1,
                       pos1=rep(LDmat$snp.id, length(LDmat$snp.id)),
                       pos2=rep(LDmat$snp.id, each=length(LDmat$snp.id)))

  LDlong <- LDlong[!is.na(cor)]

  setnames(LDlong, "pos1", "id")
  LDlong.ag <- merge(LDlong, target.variantId, by="id", allow.cartesian=T)
  setnames(LDlong.ag, "id", "pos1")

  setnames(LDlong.ag, "pos2", "id")
  LDlong.ag <- merge(LDlong.ag, target.variantId, by="id", allow.cartesian=T)
  setnames(LDlong.ag, "id", "pos2")

  LDlong.ag.ag <- LDlong.ag[pos1!=pos2,
                                list(r2.mean=mean(cor^2), r2.lci=quantile(cor^2, .025), r2.uci=quantile(cor^2, .975)),
                                list(bothSig=(sig.x==T & sig.y==T),
                                     sameChr=(chr.x==chr.y), term=term.x)]

  save(LDlong.ag.ag, LDlong.ag, file="~/LDlong.Rdata")


#### plot
  scp aob2x@rivanna.hpc.virginia.edu:~/LDlong.Rdata ~/.

  library(data.table)
  library(ggplot2)

  load("~/LDlong.Rdata")


  LDlong.sig.male <- LDlong.ag[term.x=="male"][sig.x==T & sig.y==T]
  LDlong.sig.male[,pos1o:=pos1]
  LDlong.sig.male[,pos2o:=pos2]
  LDlong.sig.male[pos2<pos1, pos1o:=pos2]
  LDlong.sig.male[pos2<pos1, pos2o:=pos1]


  LDlong.sig.male[,pos1.f:=as.numeric(as.factor(pos1o))]
  LDlong.sig.male[,pos2.f:=as.numeric(as.factor(pos2o))]

  table(LDlong.sig.male$pos1.f <= LDlong.sig.male$pos2.f)


  pl <-
  ggplot(data=LDlong.sig.male, aes(x=pos1.f, y=pos2.f, fill=cor^2)) +
  geom_tile()

  pl



  summary(t1 <- lm(cor~I(sig.x==T & sig.y==T) + I(chr.x==chr.y), LDlong.ag[pos1!=pos2]))



  ggplot(data=LDlong.ag.ag,
         aes(x=sameChr, y=r2.mean, group=interaction(bothSig, sameChr), color=bothSig)) +
  geom_point(position=position_dodge(width = .25)) +
  geom_errorbar(aes(ymin=r2.lci, ymax=r2.uci), width=.1,
                position=position_dodge(width = .25)) +
  facet_wrap(~term)
