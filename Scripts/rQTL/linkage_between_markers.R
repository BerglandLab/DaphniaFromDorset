### module load gcc/7.1.0  openmpi/3.1.4 R/4.1.1; R

### libraries
  library(data.table)
  library(foreach)
  library(data.table)
  library(SeqArray)

############
### data ###
############

### load QTL mapping output
  perm <- 0; set <- "AxC"
  load(paste("/scratch/aob2x/daphnia_hwe_sims/lmer4qtl/v3.perm", perm, ".set.", set, ".Rdata", sep=""))

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
  load(file=paste("/project/berglandlab/alan/lme4qtl_output.", set, ".long.Rdata", sep=""))
  load(file="/home/aob2x/lme4qtl_output.v3.AxC.long.Rdata")



########################
### get summary info ###
########################

  ### identify representative clone id per SC for phenotyped strains
    target.sampleId <- epp[,list(.N, sampleid=Clone[which.max(medrd)]), SCB]

  ### get positions of marker SNPs
    target.variantId <-   lmer.gwas[,list(.N, nPos=length(unique(pos))), list(id)]

################
### analysis ###
################

  LDmat <- snpgdsLDMat(genofile,
                      sample.id=target.sampleId$sampleid,
                      snp.id=target.variantId$id,
                      slide=0,
                      method=c("corr"), mat.trim=FALSE,
                      num.thread=1L, with.id=TRUE, verbose=TRUE)

  seqSetFilter(genofile, sample.id=target.sampleId$sampleid, variant.id=target.variantId$id)
