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
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")

#Load SNP file
  load("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/finalsetsnpset01_20200623.Rdata")
  seqSetFilter(genofile, variant.id=finalsetsnpset01)

# Set individual filter
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  scsub <- sc[LabGenerated==0 & Nonindependent==0 & Species=="pulex"]
  scsub$population <- str_replace(scsub$population, "Dcat", "DCat")
  scsub <- scsub[population=="DCat" | population=="D8" | population=="DBunk" | population=="D10"]
  setkey(scsub, clone)
  scsubids <- scsub$clone

  seqSetFilter(genofile, sample.id=scsubids)

### set some global parameters
  maf <- 0.001
  missing.rate <- 0.15
  threads <- 10

  ibs <- snpgdsIBS(genofile, snp.id=finalsetsnpset01, sample.id=scsubids, num.thread=20, maf=maf,
        missing.rate=0.15, autosome.only = FALSE)
  #Working space: 498 samples, 114,711 SNVs

  ### a bit of re-formating of the ibs matrix
		ibs.mat <- ibs$ibs
		rownames(ibs.mat) <- ibs$sample.id
		colnames(ibs.mat) <- ibs$sample.id

    ibs.matdt <- as.data.table(ibs.mat)
		ibs.matdt$cloneA <- scsubids
		ibs.long<- melt(ibs.matdt, measure.vars=scsubids, variable.name="cloneB", value.name="IBS")

  		### remake ibs.long
  				ibs.matdt <- as.data.table(ibs.mat)
  				ibs.matdt$cloneA <- scsubids
  				ibs.long<- melt(ibs.matdt, measure.vars=scsubids, variable.name="cloneB", value.name="IBS")
  				ibs.long <- na.omit(ibs.long)

  		### first, need to tack in SC identities
  			sctmp <- sc
  			#sctmp$SC <- ifelse(sc$clone=="April_2017_DCat_5", "B", sctmp$SC)
  			sctmp$population <- ifelse(sctmp$population=="Dcat", "DCat", sctmp$population)
  			sctmp$clone <- ifelse(sctmp$clone=="March20_2018_Dcat_1", "March20_2018_DCat_1", sctmp$clone)
  			sctmp$clone <- ifelse(sctmp$clone=="March20_2018_Dcat_2", "March20_2018_DCat_2", sctmp$clone)


  			setnames(ibs.long, "cloneA", "clone")
  			setkey(ibs.long, "clone")
  			setkey(sctmp, "clone")
  			ibs.long <- merge(ibs.long, sctmp)
  			setnames(ibs.long, "clone", "cloneA")
  			setnames(ibs.long, "SC", "SC.A")

  			setnames(ibs.long, "cloneB", "clone")
  			setkey(ibs.long, "clone")
  			setkey(sctmp, "clone")
  			ibs.long <- merge(ibs.long, sctmp)
  			setnames(ibs.long, "clone", "cloneB")
  			setnames(ibs.long, "SC", "SC.B")

  			ibs.long <- ibs.long[Species.x=="pulex" & Species.y=="pulex"]
  			ibs.long <- ibs.long[population.x!="Dramp" & population.y!="Dramp" & population.x!="DLily" &
  				population.y!="DLily" & population.x!="DOil" & population.y!="DOil" & population.x!="DMud" &
  				population.y!="DMud"]

  			### re-rank SCs based on size in pond

  			ibs.long <- ibs.long[,c("cloneA", "cloneB", "SC.A", "SC.B", "IBS"), with=F]
  							ibs.long[,pondA := tstrsplit(cloneA, "_")[[3]]]
  							ibs.long[,pondB := tstrsplit(cloneB, "_")[[3]]]


  							ibs.long.ag <- ibs.long[,list(pond.n = length(IBS)), list(pondA, SC.A) ]
  							ibs.long.ag$pond.n<- ifelse(ibs.long.ag$SC.A=="OO", 1, ibs.long.ag$pond.n)
  							setkey(ibs.long.ag, pondA, pond.n)
  							ibs.long.ag[,pond.sc.rank := ibs.long.ag[,list(pond.sc.rank = rank(-pond.n, ties="random")), list(pondA)]$pond.sc.rank]
  							ibs.long.ag[,pond.sc.rank := letters[pond.sc.rank]]

  		### be lazy and write a loop
  							ibs.long.2 <- foreach(i=1:dim(ibs.long.ag)[1], .combine="rbind")%do%{
  											temp <- ibs.long[pondA==ibs.long.ag$pondA[i] & SC.A==ibs.long.ag$SC.A[i]]
  											temp[,sc.a:=ibs.long.ag$pond.sc.rank[i]]
  											temp
  							}

  							ibs.long.3 <- foreach(i=1:dim(ibs.long.ag)[1], .combine="rbind")%do%{
  											temp <- ibs.long.2[pondB==ibs.long.ag$pondA[i] & SC.B==ibs.long.ag$SC.A[i]]
  											temp[,sc.b:=ibs.long.ag$pond.sc.rank[i]]
  											temp
  							}

  							ibs.long <- ibs.long.3

  							#setkey(ibs.long, pondA, pondB, sc.a, sc.b)

  							### next, generate [s]uper[c]lone[i]ds for individual 'A' and 'B'

  					    ibs.long[,scid.a := paste(sc.a, sprintf("%03d", as.numeric(as.factor(cloneA))), sep=".")]
  					    ibs.long[,scid.b := paste(sc.b, sprintf("%03d", as.numeric(as.factor(cloneB))), sep=".")]

  					### group on pond

  					    ibs.long[,pondA := factor(pondA, levels=c("D10", "DCat", "D8", "DBunk"))]
  					    ibs.long[,pondB := factor(pondB, levels=c("D10", "DCat", "D8", "DBunk"))]

  					    ibs.long[,scid.a := paste(LETTERS[as.numeric(pondA)], scid.a, sep=".")]
  					    ibs.long[,scid.b := paste(LETTERS[as.numeric(pondB)], scid.b, sep=".")]

  					    ibs.long[,scid.a := as.factor(scid.a)]
  					    ibs.long[,scid.b := as.factor(scid.b)]

  							### tack in buffer cloneIds for graphical purposes
  					        ibs.long <- rbind(ibs.long,
  					                          data.table(scid.a=c(paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk",
  					                                "DCat")]$pondA)], min(ibs.long$SC.A), sep=".")), c("000"), sep="."),
  					                                 paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk",
  					                                 "DCat")]$pondA)], max(ibs.long$SC.A), sep=".")), c("999"), sep=".")),
  					                                     scid.b=c(paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk",
  					                                     "DCat")]$pondB)], min(ibs.long$SC.B), sep=".")), c("000"), sep="."),
  					                                              paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk",
  					                                              "DCat")]$pondB)], max(ibs.long$SC.B), sep=".")), c("999"), sep="."))),
  					                                                fill=T)
  					         ibs.long[,scid.a := factor(scid.a, levels=sort(unique(as.character(scid.a))))]
  					         ibs.long[,scid.b := factor(scid.b, levels=sort(unique(as.character(scid.b))))]

  					    ### make lower triangle poofy-de-poof
  					        ibs.long[,dist.noTri := IBS]
  					        ibs.long[as.numeric(scid.a)>as.numeric(scid.b), dist.noTri:=NA]

  									### make pond bounding boxes
  											ibs.long.ag <- data.table(scid.a.min=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "DCat")]$pondA)], min(ibs.long$SC.A, na.rm=T), sep=".")), c("000"), sep="."),
  																								scid.a.max=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondA%in%c("D10", "D8", "DBunk", "DCat")]$pondA)], max(ibs.long$SC.A, na.rm=T), sep=".")), c("999"), sep="."),
  																								scid.b.min=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "DCat")]$pondB)], min(ibs.long$SC.B, na.rm=T), sep=".")), c("000"), sep="."),
  																								scid.b.max=paste(unique(paste(LETTERS[as.numeric(ibs.long[pondB%in%c("D10", "D8", "DBunk", "DCat")]$pondB)], max(ibs.long$SC.B, na.rm=T), sep=".")), c("999"), sep="."))

  											ibs.long.ag[,scid.a.min := as.numeric(factor(scid.a.min, levels=sort(unique(as.character(ibs.long$scid.a)))))]
  											ibs.long.ag[,scid.a.max := as.numeric(factor(scid.a.max, levels=sort(unique(as.character(ibs.long$scid.a)))))]
  											ibs.long.ag[,scid.b.min := as.numeric(factor(scid.b.min, levels=sort(unique(as.character(ibs.long$scid.b)))))]
  											ibs.long.ag[,scid.b.max := as.numeric(factor(scid.b.max, levels=sort(unique(as.character(ibs.long$scid.b)))))]


                    library(ggplot2)
                    library(viridis)

                    load("ibs.long.Rdata")

  									### plot it
  											h.just <- .25
  											v.just <- .25
  											l.size <- 1.5
  										 corrmatrix <- ggplot(data=ibs.long, aes(scid.a, scid.b, fill=IBS)) +
  											geom_raster() +
  											scale_fill_viridis(option="D")

                    ggsave(corrmatrix, file="corrmatrix_20201020.pdf")
