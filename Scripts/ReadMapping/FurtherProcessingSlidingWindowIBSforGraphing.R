#!/usr/bin/env Rscript

    ### libraries
      library(data.table)
      library(foreach)
      library(ggplot2)
      library(tidyverse)
      library(doMC)
      registerDoMC(20)


  ### Load processed IBS by window files
      inputfiles <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="mscACvsrest.ag")

      totalibsag <- foreach(i=1:length(inputfiles), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfiles[i], sep="")
        load(f)
        mscACvsrest.ag
      }

      inputfilesPO <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="mscACvsrest_POBC")

      totalibsagPO <- foreach(i=1:length(inputfilesPO), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesPO[i], sep="")
        load(f)
        mscACvsrest_POBC
      }

      totalibsagPO$ACvspond <- ifelse(totalibsagPO$ACvspond=="A_obtusa" & totalibsagPO$SCA_C=="C", "C_obtusa", ifelse(
        totalibsagPO$ACvspond=="A_pulicaria" & totalibsagPO$SCA_C=="C", "C_pulicaria", totalibsagPO$ACvspond
      ))

      inputfilesDBunkD8 <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="D8DBunkDCat.ag")

      totalibsagD8DBunk <- foreach(i=1:length(inputfilesDBunkD8), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesDBunkD8[i], sep="")
        load(f)
        D8DBunkDCat.ag
      }


      inputfilesObtusaPulicaria <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="D8DBunkDCatvsPO.ag")

      totalibsObPul <- foreach(i=1:length(inputfilesObtusaPulicaria), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesObtusaPulicaria[i], sep="")
        load(f)
        D8DBunkDCatvsPO.ag
      }

      inputfilesD10W <- list.files(path="/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", pattern="D8DBunkDCatvsD10W.ag")

      totalibsD10W <- foreach(i=1:length(inputfilesD10W), .combine="rbind")%do%{
        f=paste("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/", inputfilesD10W[i], sep="")
        load(f)
        D8DBunkDCatvsD10W.ag
      }

  ### Prep files for analysis
      total <- rbind(totalibsag, totalibsagPO, totalibsagD8DBunk, totalibsD10W, totalibsObPul)

      total$ACvspondB <- ifelse(total$ACvspond=="D8_D10" | total$ACvspond=="DBunk_D10" |
        total$ACvspond=="DCat_D10", "D8DBunkDCatvsD10", ifelse(total$ACvspond=="D8_W1" |
        total$ACvspond=="DBunk_W1" | total$ACvspond=="DCat_W1", "D8DBunkDCatvsW1", ifelse(
        total$ACvspond=="D8_W6" | total$ACvspond=="DBunk_W6" | total$ACvspond=="DCat_W6",
        "D8DBunkDCatvsW6", total$ACvspond
      )))

      total$type <- ifelse(total$ACvspondB=="A_D8DBunk_C_D8DBunk", "AvsC", ifelse(
        total$ACvspondB=="D8_D8" | total$ACvspondB=="DBunk_DBunk" | total$ACvspondB=="DCat_DCat",
        "withinpond", ifelse(total$ACvspondB=="D8_DBunk" | total$ACvspondB=="D8_DCat" |
        total$ACvspondB=="DBunk_DCat", "betweenpond", ifelse(total$ACvspondB=="A_D10" |
        total$ACvspondB=="A_W1" | total$ACvspondB=="A_W6" | total$ACvspondB=="A_pulicaria" |
        total$ACvspondB=="A_obtusa", "Avsout", ifelse(total$ACvspondB=="C_D10" |
        total$ACvspondB=="C_W1" | total$ACvspondB=="C_W6" | total$ACvspondB=="C_pulicaria" |
        total$ACvspondB=="C_obtusa", "Cvsout", "othervsout"
      )))))

      total$type <- factor(total$type, levels=c("AvsC", "withinpond", "betweenpond",
        "Avsout", "Cvsout", "othervsout"))
      total$ACvspondB <- factor(total$ACvspondB, levels=c("A_D8DBunk_C_D8DBunk",
        "D8_D8", "DBunk_DBunk", "DCat_DCat", "D8_DBunk", "D8_DCat", "DBunk_DCat", "A_D10",
        "A_W1", "A_W6", "A_pulicaria", "A_obtusa", "C_D10", "C_W1", "C_W6", "C_pulicaria",
        "C_obtusa", "D8DBunkDCatvsD10", "D8DBunkDCatvsW1", "D8DBunkDCatvsW6",
        "D8DBunkDCatvsPulicaria", "D8DBunkDCatvsObtusa"))

### Calculate IBS ratio by window and comparison for AvsC compared to all other focal pond comparisons
  AvsC <- total[ACvspondB=="A_D8DBunk_C_D8DBunk"]
  AvsCsub <- AvsC[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
  colnames(AvsCsub) <- c("window", "chr", "start", "stop", "AvsCmeanIBS")

  D8_D8 <- total[ACvspondB=="D8_D8"]
  D8_D8$SCcompare <- paste(D8_D8$SCA_C, D8_D8$SCB_C, sep="_vs_")
  D8D8SCcomp <- unique(D8_D8$SCcompare)

  AvsC_D8D8 <- foreach(i=1:length(D8D8SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8D8SCcomp[i]
    tmp <- D8_D8[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8D8", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  DBunk_DBunk <- total[ACvspondB=="DBunk_DBunk"]
  DBunk_DBunk$SCcompare <- paste(DBunk_DBunk$SCA_C, DBunk_DBunk$SCB_C, sep="_vs_")
  DBunkDBunkSCcomp <- unique(DBunk_DBunk$SCcompare)

  AvsC_DBunkDBunk <- foreach(i=1:length(DBunkDBunkSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- DBunkDBunkSCcomp[i]
    tmp <- DBunk_DBunk[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_DBunkDBunk", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  DCat_DCat <- total[ACvspondB=="DCat_DCat"]
  DCat_DCat$SCcompare <- paste(DCat_DCat$SCA_C, DCat_DCat$SCB_C, sep="_vs_")
  DCatDCatSCcomp <- unique(DCat_DCat$SCcompare)

  AvsC_DCatDCat <- foreach(i=1:length(DCatDCatSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- DCatDCatSCcomp[i]
    tmp <- DCat_DCat[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_DCatDCat", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  D8_DBunk <- total[ACvspondB=="D8_DBunk"]
  D8_DBunk$SCcompare <- paste(D8_DBunk$SCA_C, D8_DBunk$SCB_C, sep="_vs_")
  D8DBunkSCcomp <- unique(D8_DBunk$SCcompare)

  AvsC_D8DBunk <- foreach(i=1:length(D8DBunkSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DBunkSCcomp[i]
    tmp <- D8_DBunk[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8DBunk", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  D8_DCat <- total[ACvspondB=="D8_DCat"]
  D8_DCat$SCcompare <- paste(D8_DCat$SCA_C, D8_DCat$SCB_C, sep="_vs_")
  D8DCatSCcomp <- unique(D8_DCat$SCcompare)

  AvsC_D8DCat <- foreach(i=1:length(D8DCatSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DCatSCcomp[i]
    tmp <- D8_DCat[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8DCat", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  DBunk_DCat <- total[ACvspondB=="DBunk_DCat"]
  DBunk_DCat$SCcompare <- paste(DBunk_DCat$SCA_C, DBunk_DCat$SCB_C, sep="_vs_")
  DBunkDCatSCcomp <- unique(DBunk_DCat$SCcompare)

  AvsC_DBunkDCat <- foreach(i=1:length(DBunkDCatSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- DBunkDCatSCcomp[i]
    tmp <- DBunk_DCat[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_DBunkDCat", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  totalpropsimACother <- rbind(AvsC_D8D8, AvsC_DBunkDBunk, AvsC_DCatDCat, AvsC_D8DBunk, AvsC_D8DCat, AvsC_DBunkDCat)

  save(totalpropsimACother, file="totalpropsimACother_20200818.Rdata")

  D8DBunkDCatvsD10 <- total[ACvspondB=="D8DBunkDCatvsD10"]
  D8DBunkDCatvsD10$SCcompare <- paste(D8DBunkDCatvsD10$SCA_C, D8DBunkDCatvsD10$SCB_C, sep="_vs_")
  D8DBunkDCatvsD10SCcomp <- unique(D8DBunkDCatvsD10$SCcompare)

  AvsC_D8DBunkDCatvsD10 <- foreach(i=1:length(D8DBunkDCatvsD10SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DBunkDCatvsD10SCcomp[i]
    tmp <- D8DBunkDCatvsD10[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8DBunkDCatvsD10", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  D8DBunkDCatvsW1 <- total[ACvspondB=="D8DBunkDCatvsW1"]
  D8DBunkDCatvsW1$SCcompare <- paste(D8DBunkDCatvsW1$SCA_C, D8DBunkDCatvsW1$SCB_C, sep="_vs_")
  D8DBunkDCatvsW1SCcomp <- unique(D8DBunkDCatvsW1$SCcompare)

  AvsC_D8DBunkDCatvsW1 <- foreach(i=1:length(D8DBunkDCatvsW1SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DBunkDCatvsW1SCcomp[i]
    tmp <- D8DBunkDCatvsW1[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8DBunkDCatvsW1", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  D8DBunkDCatvsW6 <- total[ACvspondB=="D8DBunkDCatvsW6"]
  D8DBunkDCatvsW6$SCcompare <- paste(D8DBunkDCatvsW6$SCA_C, D8DBunkDCatvsW6$SCB_C, sep="_vs_")
  D8DBunkDCatvsW6SCcomp <- unique(D8DBunkDCatvsW6$SCcompare)

  AvsC_D8DBunkDCatvsW6 <- foreach(i=1:length(D8DBunkDCatvsW6SCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DBunkDCatvsW6SCcomp[i]
    tmp <- D8DBunkDCatvsW6[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8DBunkDCatvsW6", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  D8DBunkDCatvsPulicaria <- total[ACvspondB=="D8DBunkDCatvsPulicaria"]
  D8DBunkDCatvsPulicaria$SCcompare <- paste(D8DBunkDCatvsPulicaria$SCA_C, D8DBunkDCatvsPulicaria$SCB_C, sep="_vs_")
  D8DBunkDCatvsPulicariaSCcomp <- unique(D8DBunkDCatvsPulicaria$SCcompare)

  AvsC_D8DBunkDCatvsPulicaria <- foreach(i=1:length(D8DBunkDCatvsPulicariaSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DBunkDCatvsPulicariaSCcomp[i]
    tmp <- D8DBunkDCatvsPulicaria[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8DBunkDCatvsPulicaria", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  D8DBunkDCatvsObtusa <- total[ACvspondB=="D8DBunkDCatvsObtusa"]
  D8DBunkDCatvsObtusa$SCcompare <- paste(D8DBunkDCatvsObtusa$SCA_C, D8DBunkDCatvsObtusa$SCB_C, sep="_vs_")
  D8DBunkDCatvsObtusaSCcomp <- unique(D8DBunkDCatvsObtusa$SCcompare)

  AvsC_D8DBunkDCatvsObtusa <- foreach(i=1:length(D8DBunkDCatvsObtusaSCcomp), .combine="rbind")%do%{
    print(paste("Getting data: ", i, sep=""))
    sc <- D8DBunkDCatvsObtusaSCcomp[i]
    tmp <- D8DBunkDCatvsObtusa[SCcompare==sc]
    tmpsub <- tmp[, c("window", "chr", "start", "stop", "meanIBS"), with=TRUE]
    colnames(tmpsub) <- c("window", "chr", "start", "stop", "tmp")
    setkey(AvsCsub, window, chr, start, stop)
    setkey(tmpsub, window, chr, start, stop)
    m <- merge(AvsCsub, tmpsub)
    m$ratiodiff <- m$AvsCmeanIBS/m$tmp
    stats <- data.table(comppond="AvsC_D8DBunkDCatvsObtusa", compsc=sc, mean_ratiosim=mean(m$ratiodiff),
      sd_ratiosim=sd(m$ratiodiff))
  }

  totalpropsimDBunkD8DCatvsOut <- rbind(AvsC_D8DBunkDCatvsD10, AvsC_D8DBunkDCatvsW1, AvsC_D8DBunkDCatvsW6, AvsC_D8DBunkDCatvsPulicaria, AvsC_D8DBunkDCatvsObtusa)

  save(totalpropsimDBunkD8DCatvsOut, file="totalpropsimDBunkD8DCatvsOut_20200818.Rdata")

  ### Manipulate files
    totalgenomewide <- rbind(totalpropsimACother,totalpropsimDBunkD8DCatvsOut)

    totalgenomewidenoog <- totalgenomewide[comppond!="AvsC_A_obtusa" & comppond!="AvsC_A_pulicaria" &
      comppond!="AvsC_C_obtusa" & comppond!="AvsC_C_pulicaria" & comppond!="AvsC_D8DBunkDCatvsObtusa" &
      comppond!="AvsC_D8DBunkDCatvsPulicaria"]

    totalgenomewidenoognoAC <- totalgenomewidenoog[comppond!="AvsC_A_D10" & comppond!="AvsC_A_W1" &
      comppond!="AvsC_A_W6" & comppond!="AvsC_C_D10" & comppond!="AvsC_C_W1" & comppond!="AvsC_C_W6"]

    totalgenomewidenoognoAC$type <- ifelse(totalgenomewidenoognoAC$comppond=="AvsC_D8D8", "AvsC_D8D8", ifelse(
      totalgenomewidenoognoAC$comppond=="AvsC_DBunkDBunk" | totalgenomewidenoognoAC$comppond=="AvsC_DCatDCat",
      "AvsC_WithinPond", ifelse(totalgenomewidenoognoAC$comppond=="AvsC_D8DBunk" |
      totalgenomewidenoognoAC$comppond=="AvsC_D8DCat" | totalgenomewidenoognoAC$comppond=="AvsC_DBunkDCat",
      "AvsC_BetweenPond", totalgenomewidenoognoAC$comppond
    )))

    totalgenomewidenoognoAC$type <- factor(totalgenomewidenoognoAC$type, levels=c("AvsC_D8D8", "AvsC_WithinPond",
      "AvsC_BetweenPond", "AvsC_D8DBunkDCatvsD10", "AvsC_D8DBunkDCatvsW1", "AvsC_D8DBunkDCatvsW6"))

    totalgenomewidenoognoACsub <- totalgenomewidenoognoAC[, c("comppond", "mean_ratiosim", "sd_ratiosim", "type"), with=FALSE]

    save(totalgenomewidenoognoACsub, file="totalgenomewidenoognoACsub.Rdata")
