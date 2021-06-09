#module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(SeqArray)

### make SNP table
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile))
  setkey(snp.dt, id)



### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  f1s <- sc[AxCF1Hybrid==1]$clone


### focal locus: QTL 11
  chr.i <- "Scaffold_9199_HRSCAF_10755"
  pos.i <- 6229430
  id.i <- snp.dt[chr==chr.i][pos==pos.i]$id

  seqSetFilter(genofile, sample.id=f1s, variant.id=snp.dt[chr==chr.i][pos==pos.i]$id)
  seqGetData(genofile, "$dosage")

  ppl[variable==snp.dt[chr==chr.i][pos==pos.i]$id]
