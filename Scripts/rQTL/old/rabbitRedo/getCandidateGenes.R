### libraries
  library(data.table)
  #library(cowplot)
  library(foreach)
  library(SeqArray)

### load QTLSeqR output
  #load("/mnt/sammas_storage/bergland-lab/alan/harp_summarized_play.Rdata")
  load(file="/nv/vol186/bergland-lab/alan/peaks.Rdata")
  setnames(peaks, "CHROM", "chr")

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone=
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### which F1s?
  f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.delim")

### open GDS
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)


  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       key="chr")

### restrict to A/C and top QTL
  seqResetFilter(genofile)
  seqSetFilter(genofile,
               variant.id=snp.dt[chr=="Scaffold_9199_HRSCAF_10755"][pos>= 6229430-7500 & pos<=6229430+7500]$id,
               sample.id=sc[SC%in%c("A", "C")]$clone)

  dat <- seqGetData(genofile, "$dosage")

   dat <- data.table(dp=expand.grid(dat)$Var1,
                     clone=rep(seqGetData(genofile, "sample.id"), dim(dat)[2]),
                     variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(dat)[1]))

  setkey(dat, clone)
  setkey(sc, clone)

  dat <- merge(dat, sc)


  dat.ag <- dat[,list(freq.mu.A=round(mean(dp[SC=="A"]/2, na.rm=T), 1),
                      freq.mu.C=round(mean(dp[SC=="C"]/2, na.rm=T), 1)), list(variant.id)]


  mrk <- dat.ag[(freq.mu.A!=1 & freq.mu.C!=1) & (freq.mu.A!=0 & freq.mu.C!=0)]

  seqSetFilter(genofile,
                variant.id=mrk[freq.mu.A==.5 & freq.mu.C==.5]$variant.id)

  tmp <- seqGetData(genofile, "annotation/info/ANN")



  len1 <- tmp$length
  len2 <- tmp$data

  snp.dt1 <- data.table(len=rep(len1, times=len1),
                        ann=len2,
                        id=rep(snp.dt$variant.id, times=len1))

# Extracting data between the 1nd and 2nd | symbol
  snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
  snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]
  snp.dt1[,categor:=tstrsplit(snp.dt1$ann,"\\|")[[3]]]

  ann.ag <- snp.dt1[,list(.N), list(class, gene)]

  table(snp.dt1$class, snp.dt1$gene)
  ann.ag[class=="missense_variant"]

  table(snp.dt1$categor)

  ann.ag.ag <- ann.ag[,list(pnps=N[class=="missense_variant"]/N[class=="synonymous_variant"],
                            N=N[class=="missense_variant"] + N[class=="synonymous_variant"]), list(gene)]


  
# Collapsing additional annotations to original SNP vector length
  snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=",")),
                        list(variant.id=id)]

  snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
