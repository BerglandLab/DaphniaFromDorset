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

### Open genofile
  genofile <- seqOpen("MapJune2020_ann.seq.gds")

### Read in SNP file
  snpsvarPulex <- fread("snpsvarpulexpresentinhalf_table_20200623")

### Filter SNPs
  seqSetFilter(genofile, variant.id=snpsvarPulex$variant.ids)

### Pull out SNPs that are polymorphic within DCat, D8, DBunk at a minor allele frequency of 0.05
  sample.ids <- seqGetData(genofile, "sample.id")
  sampleidsdt <- as.data.table(sample.ids)

  temp <- unlist(strsplit(as.character(sampleidsdt$sample.ids), split="_"))
  mat <- matrix(temp, ncol=4, byrow=TRUE)
  matdat <- as.data.table(mat)
  sampleidsdt$population <- matdat$V3
  sampleidsdtD8DBunkDCat <- sampleidsdt[population=="D8" | population=="DBunk" | population=="DCat"]

	samplestouseB <- sampleidsdtD8DBunkDCat$sample.ids[sample.ids!="Spring_2017_DBunk_340" &
		sample.ids!="March20_2018_DBunk_39" & sample.ids!="Fall_2016_D10_54" &
		sample.ids!="March20_2018_D8_19" & sample.ids!="March15_2019_DBunk_MomPE1" &
		sample.ids!="March15_2019_DBunk_MomPE20" & sample.ids!="April_2017_Dbarb_11" &
		sample.ids!="March20_2018_DBunk_26" & sample.ids!="March20_2018_DBunk_37" &
		sample.ids!="March20_2018_DBunk_42" & sample.ids!="March20_2018_DBunk_10" &
		sample.ids!="March20_2018_DBunk_18" & sample.ids!="March20_2018_DBunk_21" &
		sample.ids!="March20_2018_DBunk_22" & sample.ids!="March20_2018_DBunk_23" &
		sample.ids!="March20_2018_DBunk_38" & sample.ids!="March20_2018_DBunk_40" &
		sample.ids!="March20_2018_DBunk_41" & sample.ids!="March20_2018_DBunk_43"]

	seqSetFilter(genofile, sample.id=samplestouseB)

	snps.dt <- data.table(variant.ids=seqGetData(genofile, "variant.id"),
			af=seqAlleleFreq(genofile, .progress=T))

### filter down to SNPs that are polymorphic in pulex and MAF of at least 0.05
	pulexD8DCatDBunkPoly <- snps.dt[af < 0.95 & af>0.05]

	pulexD8DCatDBunkPolyids <- pulexD8DCatDBunkPoly$variant.ids
	seqSetFilter(genofile, variant.id=pulexD8DCatDBunkPolyids)

	seqGDS2VCF(genofile, "polyDorsetMAF05.vcf.gz")
