# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(tidyverse)

### make feature counts table
  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam", "starReadsPerGene.out.tab", full.names=T)
  fn.i <- fn[1]

  fc <- foreach(fn.i=fn)%do%{
    message()
    dat <- fread(fn.i)[-c(1:4)]
    setnames(dat, c("V1", "V2"), c("GeneId", tstrsplit(fn.i, "/")%>%last%>%gsub("_starReadsPerGene.out.tab", "", .)))
    dat.df <- as.data.frame(dat[,2,with=F])
    rownames(dat.df) <- as.character(dat$GeneId)
    dat.df
  }
  fc <- do.call("cbind", fc)

  save(fc, file="~/featureCounts_STAR.Rdata")
