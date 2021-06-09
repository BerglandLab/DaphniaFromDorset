#!/usr/bin/env Rscript

### libraries
        library(gdsfmt)
        library(SNPRelate)
        library(data.table)
        library(ggplot2)
        library(foreach)
        library(lattice)
        library(tidyr)
        library(SeqArray)
        library(dplyr)
        library(tidyverse)

#Load genotype file
  genofile <- seqOpen("MapJune2020_ann.seq.gds")

#Load SNP file
  load("finalsetsnpset01_20200623.Rdata")
  seqSetFilter(genofile, variant.id=finalsetsnpset01)

# Set individual filter
  sample.ids <- seqGetData(genofile, "sample.id")
  sampleidsdt <- as.data.table(sample.ids)

### Removing lab generated clones
  temp <- unlist(strsplit(as.character(sampleidsdt$sample.ids), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  sampleidsdt$population <- matdat$V3
  sampleidsdt$population <- ifelse(sampleidsdt$population=="Dcat", "DCat", sampleidsdt$population)
  sampleidsdtsub <- sampleidsdt[population=="D8" | population=="DBunk" | population=="DCat" |
    population=="D10" | population=="DLily" | population=="DMud" | population=="DOil" |
    population=="Dramp" | population=="W1" | population=="W6" | population=="Pond21" |
    population=="Pond22" | population=="Dbarb"]

### Removing low read depth and non-independent clones
  samplestouseB <- sampleidsdtsub[sample.ids!="Spring_2017_DBunk_340" &
    sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
    sample.ids!="March20_2018_D8_19" & sample.ids!="April_2017_D8_515R" &
    sample.ids!="Lab_2019_D8_222Male" & sample.ids!="Lab_2019_D8_349Male" &
    sample.ids!="May_2017_D8_731SM" & sample.ids!="May_2017_D8_770SM" &
    sample.ids!="May_2017_D8_773SM" & sample.ids!="Spring_2017_DBunk_116SM" &
    sample.ids!="Spring_2017_DBunk_347SM" & sample.ids!="Spring_2017_DBunk_73SM" &
    sample.ids!="Spring_2016_D8_8.1"]

### Individuals listed below were not included in original superclone assignment, added in later
  samplestouseC <- samplestouseB[sample.ids!="April5_2018_D8_18105" &
    sample.ids!="April5_2018_D8_18106" & sample.ids!="April5_2018_D8_18108" &
    sample.ids!="April5_2018_D8_18109" & sample.ids!="April5_2018_D8_18111" &
    sample.ids!="April5_2018_D8_18122" & sample.ids!="March20_2018_D8_18010" &
    sample.ids!="March20_2018_D8_18025" & sample.ids!="March20_2018_D8_18028" &
    sample.ids!="March20_2018_D8_18030" & sample.ids!="March20_2018_D8_18031" &
    sample.ids!="March20_2018_DBunk_18005" & sample.ids!="March20_2018_DCat_18004" &
    sample.ids!="March20_2018_DCat_18006" & sample.ids!="March20_2018_DCat_18031" &
    sample.ids!="March_2018_DCat_18004" & sample.ids!="March20_2018_DCat_18033" &
    sample.ids!="March20_2018_DCat_18040" & sample.ids!="March20_2018_DCat_18042" &
    sample.ids!="March20_2018_DCat_18047" & sample.ids!="March20_2018_DCat_18048"]



  samplestouseCids <- samplestouseC$sample.ids

  seqSetFilter(genofile, sample.id=samplestouseCids)

### IBS

  ##IBS keeping SNPs present in 1 individuals
    ### set some global parameters
    maf <- 0.001
    missing.rate <- 0.15
    threads <- 10

    ibs <- snpgdsIBS(genofile, snp.id=finalsetsnpset01, sample.id=samplestouseC$sample.ids, num.thread=20, maf=maf,
        missing.rate=0.15, autosome.only = FALSE)

  # a bit of re-formating of the ibs matrix
    ibs.mat <- ibs$ibs
    rownames(ibs.mat) <- ibs$sample.id
    colnames(ibs.mat) <- ibs$sample.id

  # make the IBs matrix long form
    ibs.matdt <- as.data.table(ibs.mat)
    ibs.matdt$cloneA <- samplestouseC$sample.ids
    ibs.long<- melt(ibs.matdt, measure.vars=samplestouseC$sample.ids, variable.name="cloneB", value.name="IBS")

    ibs.long <- na.omit(ibs.long)

    ibs.longnoident <- ibs.long[ibs.long$cloneA!=ibs.long$cloneB]
    ibs.longnoident$cloneA <- as.factor(ibs.longnoident$cloneA)


  # Let's also remove duplicated comparisons
    ibs.longnoident$cloneAnum <- as.numeric(ibs.longnoident$cloneA)
    tmpnum <- data.table(cloneB=ibs.longnoident$cloneA, cloneBnum=ibs.longnoident$cloneAnum)
    tmpnumu <- unique(tmpnum)
    setkey(ibs.longnoident, cloneB)
    setkey(tmpnumu, cloneB)
    mibs.longnoident <- merge(ibs.longnoident, tmpnumu)
    mibs.longnoident$CloneNumComb <- ifelse(mibs.longnoident$cloneAnum > mibs.longnoident$cloneBnum,
      paste(mibs.longnoident$cloneAnum,mibs.longnoident$cloneBnum,sep="_"),
      paste(mibs.longnoident$cloneBnum,mibs.longnoident$cloneAnum,sep="_"))
    setkey(mibs.longnoident, CloneNumComb)
    ibs.longnoidentunique <- unique(mibs.longnoident, by="CloneNumComb")

  # Now get back to the original three columns
    ibs.longunique <- data.table(cloneA=ibs.longnoidentunique$cloneA,
      cloneB=ibs.longnoidentunique$cloneB, IBS=ibs.longnoidentunique$IBS)


  # Plot to determine cutoff
    ggplot(data=ibs.longunique, aes(x=IBS)) + geom_histogram(binwidth=0.001) +
      geom_vline(xintercept = 0.965, color="red") +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

### Identify superclones
  # give temporary labels to super clones based on identity of clone B
  		superclone <- ibs.long[IBS>.965]
  		superclone[,SC.sub := as.numeric(as.factor(cloneB))]

  # collapse nested superclones
  		superclone.o <- foreach(sc.i=unique(superclone$SC.sub), .combine="c")%do%{
  			paste(sort(unique(superclone[cloneB%in%superclone[SC.sub==sc.i]$cloneA]$cloneA)), collapse=";")
  			}

  	  superclone.o.uniq <- unique(superclone.o)

  	  sc.dt <- foreach(i=superclone.o.uniq, .combine="rbind")%do%{
  			data.table(clone=strsplit(i, ";")[[1]],
  			superClone.size=length(strsplit(i, ";")[[1]]),
  			superClone.index=which(i==superclone.o.uniq))
  			}

      sc.dt[,superClone.sizeRank := as.numeric(as.factor(rank(-superClone.size, ties="average")))]

  		sc.dt <- as.data.table(sc.dt %>%
  				mutate(SCnum = group_indices_(sc.dt, .dots=c("superClone.sizeRank", "superClone.index"))))

  					#temp <- unlist(strsplit(sc.dt$clone, split="_"))
  					#mat <- matrix(temp, ncol=4, byrow=TRUE)
  					#matdat <- as.data.table(mat)
  					#sc.dt$population <- matdat$V3
  					#sc.dt$year <- matdat$V2

  					save(sc.dt, file="tmp_sc.dt_20200623.Rdata")

  					### do plot to make sure our head is screwed on correctly
  					plot(superClone.sizeRank ~ superClone.size, sc.dt)

  					### label superclones with letters. What do you do when you have more than 26 superclones?
  							sc.dt[,SC:=LETTERS[SCnum]]
  							sc.dt$SC <- ifelse(sc.dt$SCnum==27, "AA", ifelse(sc.dt$SCnum==28, "AB", ifelse(sc.dt$SCnum==29, "AC",
  									ifelse(sc.dt$SCnum==30, "AD", ifelse(sc.dt$SCnum==31, "AE", ifelse(sc.dt$SCnum=="32", "AF", sc.dt$SC))))))
  							sc.dt$SC <- ifelse(sc.dt$SCnum==33, "AG", ifelse(sc.dt$SCnum==34, "AH", ifelse(sc.dt$SCnum==35, "AI",
  									ifelse(sc.dt$SCnum==36, "AJ", ifelse(sc.dt$SCnum==37, "AK", ifelse(sc.dt$SCnum=="38", "AL", sc.dt$SC))))))
  							sc.dt$SC <- ifelse(sc.dt$SCnum==39, "AM", ifelse(sc.dt$SCnum==40, "AN", ifelse(sc.dt$SCnum==41, "AO",
  									ifelse(sc.dt$SCnum==42, "AP", ifelse(sc.dt$SCnum==43, "AQ", sc.dt$SC)))))

  					### rename singleton individuals to "OO" to follow Karen's convention
  							sc.dt[superClone.size==1, SC:="OO"]

  ### This provided the base of the superclone file. 
