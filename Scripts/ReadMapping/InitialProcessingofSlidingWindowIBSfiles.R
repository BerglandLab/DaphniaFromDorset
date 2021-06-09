  #!/usr/bin/env Rscript

  ### libraries
    library(data.table)
    library(foreach)
    library(ggplot2)
    library(tidyverse)

  ############
  ### args ###
  ############

  args=commandArgs(trailingOnly=TRUE)
  varA=args[1]

   varA

  ### Let's work through one input file first
    load(paste("m_IBSbyslidingwindow_250000_10000_withpulicariaandobtusa_20200629_", varA, ".Rdata", sep=""))

  # Load superclone file
    sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
    sc$population <- str_replace(sc$population, "Dcat", "DCat")

  # Add in superclone information
    scsubA <- data.table(cloneA=sc$clone, SCA=sc$SC, popA=sc$population, yearA=sc$year, medrdA=sc$medrd, speciesA=sc$Species, sexA=sc$Sex)
    scsubB <- data.table(cloneB=sc$clone, SCB=sc$SC, popB=sc$population, yearB=sc$year, medrdB=sc$medrd, speciesB=sc$Species, sexB=sc$Sex)
    setkey(scsubA, cloneA)
    setkey(scsubB, cloneB)
    setkey(m, cloneB)
    mtmp <- merge(m, scsubB)
    setkey(mtmp, cloneA)
    msc <- merge(mtmp, scsubA)

  # Now add in new composite superclone info columns, start with making new OO SCs
    msc$SCA_B <- ifelse(msc$SCA=="OO", paste(msc$SCA, msc$cloneA, sep="_"), msc$SCA)
    msc$SCB_B <- ifelse(msc$SCB=="OO", paste(msc$SCB, msc$cloneB, sep="_"), msc$SCB)
    msc$compinfo <- paste(msc$SCA_B, msc$popA, msc$SCB_B, msc$popB, sep="_")
    msc$SCcomp <- paste(msc$SCA, msc$SCB, sep="_")

  # Now I will make two separate files, first I will make one file with all AxC comparisons, then a file that contains no A and C, and a file that contains only A and C compared to everyone else
    mscACcomp <- msc[SCcomp=="A_C" | SCcomp=="C_A"]
    mscnoAC <- msc[SCA!="A" & SCA!="C" & SCB!="A" & SCB!="C"]
    mscnoACcomp <- msc[SCcomp!="A_C" | SCcomp!="C_A"]
    mscACvsrest <- mscnoACcomp[SCA=="A" | SCA=="C" | SCB=="A" | SCB=="C"]
    mscothervsObtusaPulicaria <- mscnoAC[SCA=="G" | SCA=="P" | SCB=="G" | SCB=="P"]
    mscothervsObtusaPulicaria <- mscothervsObtusaPulicaria[SCcomp!="P_P" & SCcomp!="G_G"]

  # First aggregate the A vs C comparisons
    mscACcomp <- mscACcomp[popA!="DMud" & popB!="DMud"]
    mscACcomp$SCcomp <- ifelse(mscACcomp$SCcomp=="C_A", "A_C", mscACcomp$SCcomp)
    mscACcomp.ag <- mscACcomp[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
      list(window, chr, start, stop, speciesA, speciesB)]
    mscACcomp.ag$SCA_C=c("A")
    mscACcomp.ag$SCB_C=c("C")
    mscACcomp.ag$ACvspond=c("A_D8DBunk_C_D8DBunk")

    # Now let's aggregate the file comparing A and C to everyone else
    # Start with A versus everyone
      mscACvsrest_A <- mscACvsrest[SCA=="A" | SCB=="A"]
    # Just get's too complicated, try pulling out just A vs Ws and D10
      mscACvsrest_AW <- mscACvsrest_A[popA=="W1" | popA=="W6" | popA=="D10" | popB=="W1" | popB=="W6" | popB=="D10"]
      mscACvsrest_AW$ACvspond <- ifelse(mscACvsrest_AW$popA=="W1" | mscACvsrest_AW$popB=="W1", "A_W1", ifelse(
        mscACvsrest_AW$popA=="W6" | mscACvsrest_AW$popB=="W6", "A_W6", "A_D10"))
      mscACvsrest_AW$SCA_C <- c("A")
      mscACvsrest_AW$SCB_C <- ifelse(mscACvsrest_AW$SCB_B=="A", mscACvsrest_AW$SCA_B, mscACvsrest_AW$SCB_B)

      mscACvsrest_AW.ag <- mscACvsrest_AW[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
        list(window, chr, start, stop, speciesA, speciesB, SCA_C, SCB_C, ACvspond)]

      mscACvsrest_PO <- mscACvsrest_A[speciesA=="obtusa" | speciesB=="obtusa" | speciesA=="pulicaria" | speciesB=="pulicaria"]
      mscACvsrest_PO$ACvspond <- ifelse(mscACvsrest_PO$speciesA=="obtusa" | mscACvsrest_PO$speciesB=="obtusa", "A_obtusa", "A_pulicaria")
      mscACvsrest_PO$SCA_C <- c("A")
      mscACvsrest_PO$SCB_C <- ifelse(mscACvsrest_PO$SCB_B=="A", mscACvsrest_PO$SCA_B, mscACvsrest_PO$SCB_B)

      mscACvsrest_PO.ag <- mscACvsrest_PO[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
          list(window, chr, start, stop, SCA_C, SCB_C, ACvspond)]
      mscACvsrest_PO.ag$speciesA <- c("pulex")
      mscACvsrest_PO.ag$speciesB <- ifelse(mscACvsrest_PO.ag$ACvspond=="A_pulicaria", "pulicaria", "obtusa")

    # Do the same with C
    # Start with C versus everyone
      mscACvsrest_C <- mscACvsrest[SCA=="C" | SCB=="C"]
    # Just get's too complicated, try pulling out just A vs Ws and D10
      mscACvsrest_CW <- mscACvsrest_C[popA=="W1" | popA=="W6" | popA=="D10" | popB=="W1" | popB=="W6" | popB=="D10"]
      mscACvsrest_CW$ACvspond <- ifelse(mscACvsrest_CW$popA=="W1" | mscACvsrest_CW$popB=="W1", "C_W1", ifelse(
        mscACvsrest_CW$popA=="W6" | mscACvsrest_CW$popB=="W6", "C_W6", "C_D10"))
      mscACvsrest_CW$SCA_C <- c("C")
      mscACvsrest_CW$SCB_C <- ifelse(mscACvsrest_CW$SCB_B=="C", mscACvsrest_CW$SCA_B, mscACvsrest_CW$SCB_B)

      mscACvsrest_CW.ag <- mscACvsrest_CW[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
        list(window, chr, start, stop, speciesA, speciesB, SCA_C, SCB_C, ACvspond)]

      mscACvsrest_POB <- mscACvsrest_C[speciesA=="obtusa" | speciesB=="obtusa" | speciesA=="pulicaria" | speciesB=="pulicaria"]
      mscACvsrest_POB$ACvspond <- ifelse(mscACvsrest_POB$speciesA=="obtusa" | mscACvsrest_POB$speciesB=="obtusa", "A_obtusa", "A_pulicaria")
      mscACvsrest_POB$SCA_C <- c("C")
      mscACvsrest_POB$SCB_C <- ifelse(mscACvsrest_POB$SCB_B=="C", mscACvsrest_POB$SCA_B, mscACvsrest_POB$SCB_B)

      mscACvsrest_POB.ag <- mscACvsrest_POB[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
          list(window, chr, start, stop, SCA_C, SCB_C, ACvspond)]
      mscACvsrest_POB.ag$speciesA <- c("pulex")
      mscACvsrest_POB.ag$speciesB <- ifelse(mscACvsrest_POB.ag$ACvspond=="A_pulicaria", "pulicaria", "obtusa")

      mscACvsrest.ag <- rbind(mscACvsrest_AW.ag, mscACvsrest_CW.ag, mscACcomp.ag)

      mscACvsrest_POBC <- rbind(mscACvsrest_PO.ag, mscACvsrest_POB.ag)

      save(mscACvsrest.ag, file=paste("mscACvsrest.ag_", varA, ".Rdata", sep=""))

      save(mscACvsrest_POBC, file=paste("mscACvsrest_POBC_", varA, ".Rdata", sep=""))

### Make a file to compare all non A and C lineages from the three focal ponds to the outgroup species

          mscnoAC <- msc[SCA!="A" & SCA!="C" & SCB!="A" & SCB!="C"]
          mscothervsObtusaPulicaria <- mscnoAC[SCA=="G" | SCA=="P" | SCB=="G" | SCB=="P"]
          mscothervsObtusaPulicaria <- mscothervsObtusaPulicaria[SCcomp!="P_P" & SCcomp!="G_G"]
          mscothervsObtusaPulicaria <- mscothervsObtusaPulicaria[popA!="D10" & popB!="D10" &
            popA!="DLily" & popB!="DLily" & popA!="DOil" & popB!="DOil" & popA!="Dramp" &
            popB!="Dramp" & popA!="W1" & popB!="W1" & popA!="W6" & popB!="W6"]
          mscothervsObtusaPulicaria$ACvspond <- ifelse(mscothervsObtusaPulicaria$speciesA=="obtusa" |
            mscothervsObtusaPulicaria$speciesB=="obtusa", "D8DBunkDCatvsObtusa", "D8DBunkDCatvsPulicaria")
          mscothervsObtusaPulicaria$SCcompB <- ifelse(mscothervsObtusaPulicaria$SCA_B >
            mscothervsObtusaPulicaria$SCB_B, paste(mscothervsObtusaPulicaria$SCA_B,
            mscothervsObtusaPulicaria$SCB_B, sep="_"), paste(mscothervsObtusaPulicaria$SCB_B, mscothervsObtusaPulicaria$SCA_B, sep="_"))
          mscothervsObtusaPulicaria$SCA_C <- mscothervsObtusaPulicaria$SCA_B
          mscothervsObtusaPulicaria$SCB_C <- mscothervsObtusaPulicaria$SCB_B

          D8DBunkDCatvsPO.ag <- mscothervsObtusaPulicaria[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
            list(window, chr, start, stop, speciesA, speciesB, SCA_C, SCB_C, ACvspond)]

          save(D8DBunkDCatvsPO.ag, file=paste("D8DBunkDCatvsPO.ag_", varA, ".Rdata", sep=""))

### Pull out D10, W ponds versus three focal ponds
          mscnoACD10W <- mscnoAC[popA=="D10" | popB=="D10" | popA=="W1" | popB=="W1" | popA=="W6" | popB=="W6"]
          mscnoACD10W <- mscnoACD10W[SCA!="P" & SCB!="P" & SCA!="G" & SCB!="G"]
          mscnoACD10W <- mscnoACD10W[popB!="DLily" & popB!="DOil" & popB!="Dramp"]
          mscnoACD10W <- mscnoACD10W[SCA_B!=SCB_B]
          mscnoACD10W$SCcompB <- ifelse(mscnoACD10W$SCA_B > mscnoACD10W$SCB_B, paste(mscnoACD10W$SCA_B,
          mscnoACD10W$SCB_B, sep="_"), paste(mscnoACD10W$SCB_B, mscnoACD10W$SCA_B, sep="_"))
          mscnoACD10W <- mscnoACD10W[popA!=popB]
          mscnoACD10W$SCA_C <- ifelse(mscnoACD10W$SCA_B > mscnoACD10W$SCB_B, mscnoACD10W$SCA_B, mscnoACD10W$SCB_B)
          mscnoACD10W$SCB_C <- ifelse(mscnoACD10W$SCA_B > mscnoACD10W$SCB_B, mscnoACD10W$SCB_B, mscnoACD10W$SCA_B)
          mscnoACD10W$ACvspond <- paste(mscnoACD10W$popA, mscnoACD10W$popB, sep="_")
          mscnoACD10W$ACvspond <- ifelse(mscnoACD10W$ACvspond=="D10_D8", "D8_D10", ifelse(mscnoACD10W$ACvspond=="D10_DBunk",
          "DBunk_D10", ifelse(mscnoACD10W$ACvspond=="D10_DCat", "DCat_D10", ifelse(mscnoACD10W$ACvspond=="W1_D8", "D8_W1",
          ifelse(mscnoACD10W$ACvspond=="W6_D8", "D8_W6", ifelse(mscnoACD10W$ACvspond=="W1_DCat", "DCat_W1", ifelse(
          mscnoACD10W$ACvspond=="W6_DCat", "DCat_W6", ifelse(mscnoACD10W$ACvspond=="W1_DBunk", "DBunk_W1", ifelse(
          mscnoACD10W$ACvspond=="W6_DBunk", "DBunk_W6", mscnoACD10W$ACvspond)))))))))
          mscnoACD10W <- mscnoACD10W[ACvspond!="W1_D10" & ACvspond!="W6_D10" & ACvspond!="W6_W1"]

          D8DBunkDCatvsD10W.ag <- mscnoACD10W[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
            list(window, chr, start, stop, speciesA, speciesB, SCA_C, SCB_C, ACvspond)]

          save(D8DBunkDCatvsD10W.ag, file=paste("D8DBunkDCatvsD10W.ag_", varA, ".Rdata", sep=""))


#### Pull out just D8 and DBunk and DCat
          mscnoAC$popA <- ifelse(mscnoAC$popA=="Dcat", "DCat", mscnoAC$popA)
          mscnoAC$popB <- ifelse(mscnoAC$popB=="Dcat", "DCat", mscnoAC$popB)
          mscD8DBunk <- mscnoAC[popA=="D8" | popA=="DBunk" | popA=="DCat"]
          mscD8DBunkB <- mscD8DBunk[popB=="D8" | popB=="DBunk" | popB=="DCat"]
          mscD8DBunkB <- mscD8DBunkB[SCA_B!=SCB_B]
          mscD8DBunkB <- mscD8DBunkB[speciesA!="obtusa" & speciesB!="obtusa"]
          mscD8DBunkB$SCcompB <- ifelse(mscD8DBunkB$SCA_B > mscD8DBunkB$SCB_B, paste(mscD8DBunkB$SCA_B,
          mscD8DBunkB$SCB_B, sep="_"), paste(mscD8DBunkB$SCB_B, mscD8DBunkB$SCA_B, sep="_"))
          mscD8DBunkB$SCA_C <- ifelse(mscD8DBunkB$SCA_B > mscD8DBunkB$SCB_B, mscD8DBunkB$SCA_B, mscD8DBunkB$SCB_B)
          mscD8DBunkB$SCB_C <- ifelse(mscD8DBunkB$SCA_B > mscD8DBunkB$SCB_B, mscD8DBunkB$SCB_B, mscD8DBunkB$SCA_B)
          mscD8DBunkB$ACvspond <- paste(mscD8DBunkB$popA, mscD8DBunkB$popB, sep="_")
          mscD8DBunkB$ACvspond <- ifelse(mscD8DBunkB$ACvspond=="DBunk_D8", "D8_DBunk", ifelse(
          mscD8DBunkB$ACvspond=="DCat_D8", "D8_DCat", ifelse(mscD8DBunkB$ACvspond=="DCat_DBunk",
          "DBunk_DCat", mscD8DBunkB$ACvspond)))

          D8DBunkDCat.ag <- mscD8DBunkB[,list(meanIBS=mean(IBS), sd_IBS=sd(IBS)),
            list(window, chr, start, stop, speciesA, speciesB, SCA_C, SCB_C, ACvspond)]

          save(D8DBunkDCat.ag, file=paste("D8DBunkDCat.ag_", varA, ".Rdata", sep=""))
