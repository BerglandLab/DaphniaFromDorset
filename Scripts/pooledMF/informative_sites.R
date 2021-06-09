#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

#### makes rabbit input
#  args = commandArgs(trailingOnly=TRUE)
#  chr.i <- as.character(args[1])
#  maxcM <- as.numeric(args[2])
#  f1s.set <- as.character(args[3])
#  #chr.i <- "Scaffold_1863_HRSCAF_2081"; maxcM=10; f1s.set <- "all"

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### open GDS
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

### load in filter file
  snpFilter <- fread("snpsvarpulexpresentinhalf_table_20200623")

### make snp.dt
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       key="chr")
  setkey(snpFilter, chr, pos)
  setkey(snp.dt, chr, pos)

  snp.dt <- merge(snpFilter, snp.dt)



### make majority rule (consensus) genotype calls for
  ac.fd <- foreach(sc.i=c("A", "C"), .combine="cbind")%do%{
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sc[SC==sc.i]$clone, variant.id=snp.dt$id)

    data.table(af=seqAlleleFreq(genofile, ref.allele=1L)) ### alternate allele
  }
  setnames(ac.fd, c(1,2), c("af.A", "af.C"))
  ac.fd <- cbind(ac.fd, snp.dt)


  ac.fd[!is.na(af.A),A.geno := unlist(sapply(ac.fd[!is.na(af.A)]$af.A, function(x) c("11","12","22")[which.min(abs(x-c(0,.5,1)))]))]
  ac.fd[!is.na(af.C),C.geno := unlist(sapply(ac.fd[!is.na(af.C)]$af.C, function(x) c("11","12","22")[which.min(abs(x-c(0,.5,1)))]))]

  ac.fd[!is.na(af.A),A.delta := unlist(sapply(ac.fd[!is.na(af.A)]$af.A, function(x) min(abs(x-c(0,.5,1)))))]
  ac.fd[!is.na(af.C),C.delta := unlist(sapply(ac.fd[!is.na(af.C)]$af.C, function(x) min(abs(x-c(0,.5,1)))))]


  ac.inform <- ac.fd[(A.geno=="12" & C.geno=="11") |
                     (A.geno=="12" & C.geno=="22") |
                     (A.geno=="11" & C.geno=="12") |
                     (A.geno=="22" & C.geno=="12") |
                     (A.geno=="12" & C.geno=="12") ]

  ac.inform <- ac.inform[A.delta < 0.05 & C.delta < 0.05]
  save(ac.inform, file="/scratch/aob2x/daphnia_hwe_sims/ac_inform.Rdata")
  
### export VCF
  seqSetFilter(genofile, sample.id=sc[SC%in%c("A", "C")]$clone, variant.id=ac.inform$id)

  seqGDS2VCF(genofile,
             vcf.fn="/scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf", info.var=character(0), fmt.var=character(0))







 prop.table(table(A=ac.fd[grepl("9200", chr)][pos>8e6]$A.geno,
 C=ac.fd[grepl("9200", chr)][pos>8e6]$C.geno))


prop.table(table(ac.fd[grepl("9200", chr)][pos<8e6]$A.geno,
ac.fd[grepl("9200", chr)][pos<8e6]$C.geno))
