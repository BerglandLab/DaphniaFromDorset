
### library
  library(data.table)
  library(foreach)
  library(SeqArray)
  library(tidyverse)
  library(Rsubread)

### load GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

### ase sample table
  sampleTable <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable")
  setnames(sampleTable, "SampleName", "samp")
  sampleTable[,clone:=toupper(clone)]

### superclone file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  sc <- sc[Nonindependent==0 ]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]
  sc[,pond:=toupper(population)]

### replace clone name
  for(i in 1:8) {
    sampleTable[i,clone:=sc[grepl(sampleTable[i]$clone, clone)]$clone]
  }
  sampleTable


### snp.table
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       ref=seqGetData(genofile, "$ref"),
                       alt=seqGetData(genofile, "$alt"))
  setkey(snp.dt, chr, pos)

  ### load in filter file
    snpFilter <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/snpsvarpulexpresentinhalf_table_20200623")
    setkey(snpFilter, chr, pos)

  ## filter
    snp.dt <- merge(snpFilter, snp.dt)

### get read depth and missing data counts
  seqSetFilter(genofile, sample.id=unique(sampleTable$clone), variant.id=snp.dt$id)


  dp <- seqGetData(genofile, "annotation/format/DP")
  dosage <- seqGetData(genofile, "$dosage")

  dp.dt <- as.data.table(t(dp$data))
  setnames(dp.dt, seqGetData(genofile, "sample.id"))
  dp.dt[,id:=seqGetData(genofile, "variant.id")]
  dp.dt.l <- melt(dp.dt, id.vars="id", variable.name = "sample.id", value.name = "dp")
  dp.dt.l <- dp.dt.l[,list(dp=dp, dp.norm=dp/mean(dp, na.rm=T), id), list(sample.id)]

  dosage.dt <- as.data.table(t(dosage))
  setnames(dosage.dt, seqGetData(genofile, "sample.id"))
  dosage.dt[,id:=seqGetData(genofile, "variant.id")]
  dosage.dt.l <- melt(dosage.dt, id.vars="id", variable.name = "sample.id", value.name = "dosage")
  dosage.dt.l <- dosage.dt.l[,list(dosage=dosage, dosage.norm=dosage/mean(dosage, na.rm=T), id), list(sample.id)]

  setkey(dp.dt.l, sample.id, id)
  setkey(dosage.dt.l, sample.id, id)

  dd <- merge(dp.dt.l, dosage.dt.l)
  setkey(dd, id)
  setkey(snp.dt, id)

  dd <- merge(dd, snp.dt)
  dd[,start:=pos]
  dd[,end:=pos]

### get basic gene info
  saf <- flattenGTF("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf",
      method = "merge")

  saf <- as.data.table(saf)
  saf.dt <- saf[,list(start=min(Start),end=max(End), chr=unique(Chr)), list(GeneID)]

### foverlaps
  setkey(dd, chr, start, end)
  setkey(saf.dt, chr, start, end)
  genes <- foverlaps(saf.dt, dd)

  genes.ag <- genes[,list(dp.norm.mu=mean(dp.norm, na.rm=T),
                          missing.rate=mean(is.na(dosage)), .N),
                     list(GeneID, clone=sample.id)]

  genes.ag <- merge(genes.ag[!is.na(clone)], sampleTable, by="clone", allow.cartesian=T)

  save(genes.ag, file="~/genes_rd_missing.Rdata")
