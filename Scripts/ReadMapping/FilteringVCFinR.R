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
		library(tidyverse)

### format data properly to work with SNPRelate
        ### subsetted data for PCA stuff
                vcf.fn <- "MapJune2020_ann.vcf"
                snpgdsVCF2GDS(vcf.fn, "MapJune2020_ann.gds",
                                        method=c("biallelic.only"), snpfirstdim = FALSE)

                seqVCF2GDS(vcf.fn, "MapJune2020_ann.seq.gds")

### Open VCF and take an initial look at the SNPs
    genofile <- seqOpen("MapJune2020_ann.seq.gds")

    snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
        chr = seqGetData(genofile, "chromosome"),
        pos = seqGetData(genofile, "position"),
        dp = seqGetData(genofile, "annotation/info/DP"))

### Filter out SNPs that occur in areas flagged as having too high or low read depth when mapping the D84A 10X Illumina short reads to the D84A reference genome.
### Will also filter out SNPs on the edges of runs of Ns and at the ends of scaffolds.

    ## Read in bed file that contains regions of high/low read depth and end of stretches of Ns and ends of chromosomes
        NsChrRD <- fread("NsandDepthandChrEnd.sorted.500merged.bed")

    ## Filter out SNPs in these regions
        colnames(NsChrRD) <- c("chr", "start", "stop")
        NsChrRD$count <- c(1:26531)

        setkey(snps, chr, pos)
        initialsnps <- snps$variant.ids

        NsChrRDsnps <- foreach(i=NsChrRD$count, .combine="c")%do%{
            c=NsChrRD$chr[[i]]
            s=NsChrRD$start[[i]]
            p=NsChrRD$stop[[i]]
            temp <- snps[J(data.table(chr=c, pos=c(s:p), key="chr,pos")), nomatch=0]
            temp$variant.ids
          }

        goodsnpsnotinNsChrRD <- setdiff(initialsnps, NsChrRDsnps)

        seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRD)

        goodsnpsnotinNsChrRDtable <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
	         chr = seqGetData(genofile, "chromosome"),
	         pos = seqGetData(genofile, "position"),
	         dp = seqGetData(genofile, "annotation/info/DP"))

### Now remove SNPs in repeat masker identified regions
    ## Read in bed file that contains regions flagged by Repeat Masker
        RMout <- fread("RMoutHiCGMgoodscaff.bed")

    ## Now remove SNPs in those regions
        colnames(RMout) <- c("chr", "start", "stop")
        RMout$count <- c(1:125881)

        setkey(goodsnpsnotinNsChrRDtable, chr, pos)
        initialsnps <- goodsnpsnotinNsChrRDtable$variant.ids

        RMoutSNPs <- foreach(i=RMout$count, .combine="c")%do%{
            c=RMout$chr[[i]]
            s=RMout$start[[i]]
            p=RMout$stop[[i]]
            temp <- goodsnpsnotinNsChrRDtable[J(data.table(chr=c, pos=c(s:p), key="chr,pos")), nomatch=0]
            temp$variant.ids
          }

        goodsnpsnotinNsChrRDorRM <- setdiff(goodsnpsnotinNsChrRD, RMoutSNPs)

        seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRDorRM)

        goodsnpsnotinNsChrRDorRMtable <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
            chr = seqGetData(genofile, "chromosome"),
            pos = seqGetData(genofile, "position"),
            dp = seqGetData(genofile, "annotation/info/DP"))

### Now remove triallelic snps
    goodsnpsnotinNsChrRDorRM <- goodsnpsnotinNsChrRDorRMtable$variant.ids

    seqSetFilter(genofile, variant.id=goodsnpsnotinNsChrRDorRM)

    tri <- (seqGetData(genofile, "$num_allele"))
    tri <- as.data.table(tri)
    tri$variant.ids <- seqGetData(genofile, "variant.id")
    tri$diallelic <- ifelse(tri$tri=="2", 1, 0)

    updatesnpstouse <- tri$variant.ids[tri$tri=="2"]

    seqSetFilter(genofile, variant.id=updatesnpstouse)

    snps <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
        chr = seqGetData(genofile, "chromosome"),
        pos = seqGetData(genofile, "position"),
        dp = seqGetData(genofile, "annotation/info/DP"))

### Overall read depth filtering on snps
    quantile(snps$dp,probs=c(0.05,0.95))

  	#5%   95%
  	#3497 12115

  	dpfiltsnps <- snps[dp > 3496 & dp < 12114]
  	dpfiltsnpsids <- dpfiltsnps$variant.ids

  	seqSetFilter(genofile, variant.id=dpfiltsnpsids)

  	dpfiltsnpsdt <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
  			chr = seqGetData(genofile, "chromosome"),
  			pos = seqGetData(genofile, "position"),
  			dp = seqGetData(genofile, "annotation/info/DP"))

    write.csv(dpfiltsnpsdt, file="dpfiltsnpsdt.csv", row.names=FALSE, quote=FALSE)

### This dpfilt SNP set was used for analyses that include the outgroup species.
### Now filter to only keep SNPs fixed within Daphnia pulex
    seqSetFilter(genofile, variant.id=dpfiltsnpsids)

    ## Remove D. pulicaria, D. obtusa, and very low read depth individuals (poor quality)

    	sample.ids <- seqGetData(genofile, "sample.id")
    	sampleidsdt <- as.data.table(sample.ids)
    	samplestouseB <- sampleidsdt$sample.ids[sample.ids!="Spring_2017_DBunk_340" &
    	   sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
    	   sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
    	   sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
    	   sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
    	   sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
    	   sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
    	   sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
    	   sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
    	   sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43" &
    	   sample.ids!="2018_Pulicaria_Pond21_22" & sample.ids!="2018_Pulicaria_Pond22_21" &
    	   sample.ids!="2018_Pulicaria_Pond22_53" & sample.ids!="2018_Pulicaria_Pond22_62" &
    	   sample.ids!="2018_Pulicaria_Pond22_72"]

    	seqSetFilter(genofile, sample.id=samplestouseB)

    ## Filter down to SNPs that are polymorphic in D. pulex.
    	snps.dt <- data.table(variant.ids=seqGetData(genofile, "variant.id"),
    												af=seqAlleleFreq(genofile, .progress=T))

    	pulexPoly <- snps.dt[af!=1 & af!=0]

    	setkey(dpfiltsnps, variant.ids)
    	setkey(pulexPoly, variant.ids)

    	snpsvarPulex <- merge(dpfiltsnps, pulexPoly)

    ## Check number of individuals per SNP.
      snpstousevarPulex <- snpsvarPulex$variant.ids
  		seqSetFilter(genofile, variant.id=snpstousevarPulex)

      #Pull out genotypes

    			seqSetFilter(genofile, sample.id=samplestouseB)

    			het <- t(seqGetData(genofile, "$dosage"))
    			het <- as.data.table(het)

    			colnames(het) <- c(seqGetData(genofile, "sample.id"))
    			het$variant.ids <- seqGetData(genofile, "variant.id")

    			setkey(het, variant.ids)
    			setkey(snpsvarPulex, variant.ids)

    			mhet <- merge(snpsvarPulex, het)

    	# Set all genotypes to 1

    			mhet[mhet == 0] <- 1
    			mhet[mhet == 2] <- 1

    			mhetlong <- melt(mhet, measure.vars=samplestouseB, variable.name="clone", value.name="dosage")

    			mhetlong.ag <- mhetlong[,list(numgeno = sum(dosage, na.rm=TRUE)), list(variant.ids) ]

    			indpersnp <- ggplot(data=mhetlong.ag, aes(x=numgeno)) + geom_histogram(binwidth=10)
    			ggsave(indpersnp, file="indpersnp_20200623.pdf")

    	# Let's drop SNPs that are genotyped in less than half the individuals
    			snpsvarpulexpresentinhalf <- mhetlong.ag$variant.ids[mhetlong.ag$numgeno>283]

    			seqSetFilter(genofile, variant.id=snpsvarpulexpresentinhalf)

    			snpsvarPulex <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
    						chr = seqGetData(genofile, "chromosome"),
    						pos = seqGetData(genofile, "position"),
    						dp = seqGetData(genofile, "annotation/info/DP"))

    			write.table(snpsvarPulex, file="snpsvarpulexpresentinhalf_table_20200623", sep="\t", row.names=FALSE, quote=FALSE)

      # This SNP set represents good quality SNPs within D. pulex that are used in multiple analyses.

### Next SNP pruning. Setting maf pretty low, so SNPs are kept even if present in 1 individual.

    seqSetFilter(genofile, variant.id=snpsvarpulexpresentinhalf)

    ## Set individual filter to ensure all outgroup species are removed, low read depth and non-independent D. pulex samples are removed.
      sample.ids <- seqGetData(genofile, "sample.id")
      sampleidsdt <- as.data.table(sample.ids)

      temp <- unlist(strsplit(as.character(sampleidsdt$sample.ids), split="_"))
      mat <- matrix(temp, ncol=4, byrow=TRUE)
      matdat <- as.data.table(mat)
      sampleidsdt$population <- matdat$V3
      sampleidsdt$population <- ifelse(sampleidsdt$population=="Dcat", "DCat", sampleidsdt$population)
      sampleidsdtsub <- sampleidsdt[population=="D8" | population=="DBunk" | population=="DCat" |
        population=="D10" | population=="DLily" | population=="DMud" | population=="DOil" |
        population=="Dramp" | population=="W1" | population=="W6"]

      samplestouseB <- sampleidsdtsub[sample.ids!="Spring_2017_DBunk_340" &
      sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
      sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
      sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
      sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
      sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
      sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
      sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
      sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
      sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43" &
      sample.ids!="2018_Pulicaria_Pond21_22" & sample.ids!="2018_Pulicaria_Pond22_21" &
      sample.ids!="2018_Pulicaria_Pond22_53" & sample.ids!="2018_Pulicaria_Pond22_62" &
      sample.ids!="2018_Pulicaria_Pond22_72" & sample.ids!="April_2017_D8_515R" &
      sample.ids!="Lab_2019_D8_222Male" & sample.ids!="Lab_2019_D8_349Male" &
      sample.ids!="May_2017_D8_731SM" & sample.ids!="May_2017_D8_770SM" &
      sample.ids!="May_2017_D8_773SM" & sample.ids!="Spring_2017_DBunk_116SM" &
      sample.ids!="Spring_2017_DBunk_347SM" & sample.ids!="Spring_2017_DBunk_73SM" &
      sample.ids!="Spring_2016_D8_8.1" & sample.ids!="March_2018_D8_18030" &
      sample.ids!="March_2018_DCat_18004"]

      samplestouseBids <- samplestouseB$sample.ids

      seqSetFilter(genofile, sample.id=samplestouseBids)

  ## set some global parameters
      maf <- 0.001
      missing.rate <- 0.15
      threads <- 10

  ## SNP pruning (removing lab generated clones)
      set.seed(10000)
      snpset01 <- snpgdsLDpruning(genofile, snp.id=snpsvarpulexpresentinhalf, sample.id=samplestouseB,
      	autosome.only=FALSE, maf=maf,missing.rate=missing.rate, slide.max.bp=500, ld.threshold=0.1)
      finalsetsnpset01 <-unlist(snpset01[c(1:62)])

      seqSetFilter(genofile, variant.id=finalsetsnpset01)

      snpsvarPulexLDprune <- data.table(variant.ids = seqGetData(genofile, "variant.id"),
          chr = seqGetData(genofile, "chromosome"),
          pos = seqGetData(genofile, "position"),
          dp = seqGetData(genofile, "annotation/info/DP"))

      write.table(snpsvarPulexLDprune, file="finalsetsnpset01pulex_table_20200623", sep="\t", row.names=FALSE, quote=FALSE)

  ## This SNP set is the linkage pruned SNP set of SNPs polymorphic within D. pulex, which are used for some analyses.
