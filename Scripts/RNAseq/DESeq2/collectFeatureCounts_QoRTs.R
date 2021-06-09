# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(tidyverse)

### make feature counts table
  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam", "starReadsPerGene.out.tab", full.names=T)
  fn <- list.files("/Users/alanbergland/qorts/", full.names=T)
  fn.i <- fn[1]

  fc <- foreach(fn.i=fn)%do%{
    message()
    dat <- fread(paste(fn.i, "/QC.geneCounts.formatted.for.DESeq.txt.gz", sep=""))

    setnames(dat, c("V1", "V2"), c("GeneId", tstrsplit(fn.i, "/")%>%last))
    dat <- dat[grepl("Daphnia", GeneId)]
    dat.df <- as.data.frame(dat[,2,with=F])
    rownames(dat.df) <- as.character(dat$GeneId)
    dat.df
  }
  fc <- do.call("cbind", fc)

  save(fc, file="~/featureCounts_QoRTs.Rdata")
