#!/usr/bin/env Rscript

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(20)
  library(ggplot2)
  library(SNPRelate)
  #library(tidyverse)

  ############
  ### args ###
  ############

  args=commandArgs(trailingOnly=TRUE)
  varA=args[1]
  varB=args[2]
  varC=args[3]

### open genotype file
genofile <- seqOpen("/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann.seq.gds")

###Use non LD pruned dataset
  load("dpfiltsnps_20200623.Rdata")
  filtsnptb <- dpfiltsnps
  colnames(filtsnptb) <- c("oldvariantids", "chr", "pos", "olddp")

### Filter SNPs
  snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
      chr = seqGetData(genofile, "chromosome"),
      pos = seqGetData(genofile, "position"))
  setkey(snps, chr, pos)
  setkey(filtsnptb, chr, pos)
  msnps <- merge(filtsnptb, snps)
  msnps[,final.use:=T]

### Load superclone file
  sc <- fread("/scratch/kbb7sh/Daphnia/MappingDecember2019/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### Do some filtering on read depth and independence
  scrd5 <- sc[medrd>4 & Nonindependent==0 & LabGenerated==0 & clone!="April5_2018_D8_Male2"]
  subscrd5 <- scrd5[, c("clone", "SC", "population", "year", "Sex", "Species", "medrd"), with=FALSE]
  clonestouse <- subscrd5
  clonestouseids <- clonestouse$clone

### make windows
  setkey(msnps, chr, pos)

  chr.ag <- msnps[,list(start = min(pos), stop=max(pos)), list(chr)]
  chr.ag$length <- chr.ag$stop-chr.ag$start
  chrtouse <- chr.ag$chr[chr.ag$length>1000000]

  window.size <- 250000
  step.size <- 10000
  wins <- foreach(chr.i=chrtouse, .combine="rbind")%dopar%{
    #chr.i <- "Scaffold_1863_HRSCAF_2081"
    snp.dt.ag <- msnps[J(chr.i)][(final.use), list(start=min(pos), stop=max(pos))]

    data.table(chr=chr.i, start=seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size),
                          stop =seq(from=snp.dt.ag$start,
                                    to=snp.dt.ag$stop - window.size,
                                    by=step.size) + window.size)
  }

  ### Maybe look at how many SNPs go into each window/tree? Perhaps there are some low ones that are throwing things off?
  msnpcount <- foreach(win.i = c(1:dim(wins)[1]), .combine="rbind", .errorhandling="remove")%dopar%{

    print(paste("Getting data: ", win.i, sep=""))

    snp.dt.tmp <- msnps[J(wins[win.i]$chr)][pos>=wins[win.i]$start & pos<=wins[win.i]$stop][(final.use)]

    tmp <- data.table(window=win.i, numsnps=dim(snp.dt.tmp)[1])
    tmp

  }

  wins$numsnps <- msnpcount$numsnps
  winstouse <- wins[numsnps>1000]


### run windows
  #m <- foreach(win.i = c(1:10), .combine="rbind", .errorhandling="remove")%dopar%{
  m <- foreach(win.i =varA:varB, .errorhandling="remove", .combine="rbind")%dopar%{

    print(paste("Getting data: ", win.i, sep=""))

    snp.dt.tmp <- msnps[J(winstouse[win.i]$chr)][pos>=winstouse[win.i]$start & pos<=winstouse[win.i]$stop][(final.use)]

    seqResetFilter(genofile)
    seqSetFilter(genofile, variant.id=snp.dt.tmp$variant.ids, sample.id=clonestouseids)

    #First lets try IBS with we keep SNPs present in 1 individuals
		### set some global parameters
				maf <- 0.001
				missing.rate <- 0.15
				threads <- 10

				ibs <- snpgdsIBS(genofile, snp.id=snp.dt.tmp$variant.ids, sample.id=clonestouseids, num.thread=20, maf=maf,
							missing.rate=0.15, autosome.only = FALSE)

			### a bit of re-formating of the ibs matrix
				ibs.mat <- ibs$ibs
				rownames(ibs.mat) <- ibs$sample.id
				colnames(ibs.mat) <- ibs$sample.id

			### make the IBs matrix long form
        ibs.matdt <- as.data.table(ibs.mat)
        setkey(clonestouse, clone)
        ibs.matdt$cloneA <- clonestouse$clone
        ibs.long<- melt(ibs.matdt, measure.vars=clonestouseids, variable.name="cloneB", value.name="IBS")
				ibs.long <- na.omit(ibs.long)

        # First let's remove all identical comparisons from ibs.long
  				ibs.longnoident <- ibs.long[ibs.long$cloneA!=ibs.long$cloneB]
          ibs.longnoident$cloneA <- as.factor(ibs.longnoident$cloneA)

  			#Let's also remove duplicated comparisons - how to do this?
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
  			#Now get back to the original three columns
  				ibs.longunique <- data.table(cloneA=ibs.longnoidentunique$cloneA,
  					cloneB=ibs.longnoidentunique$cloneB, IBS=ibs.longnoidentunique$IBS)

        ibs.longunique$chr <- c(winstouse[win.i]$chr)
        ibs.longunique$start <- c(winstouse[win.i]$start)
        ibs.longunique$stop <- c(winstouse[win.i]$stop)
        ibs.longunique$window <- c(win.i)
        ibs.longunique

  }

   save(m, file=paste("m_IBSbyslidingwindow_250000_10000_withpulicariaandobtusa_20200629_", varC, ".Rdata", sep=""))

