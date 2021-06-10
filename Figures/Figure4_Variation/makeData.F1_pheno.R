#### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)

### This is the same procesing of phenotype data that goes into the lme4qtl analysis.
### copied from `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/lme4qtl.gwas.R`

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

### save F1 data object
  save(epp, male, file="/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanFigures/Fig_4_polygenic/F1_pheno.Rdata")
