#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

#scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/popPhase/trees/region.fasta ~/.


  #args = commandArgs(trailingOnly=TRUE)
  #step.bp=as.numeric(args[[1]])
  #window.bp=as.numeric(args[[2]])
#
  step.bp <- 50000
  window.bp <- 250000

### libraries
  #library(ape)
  library(data.table)
  library(foreach)

### chrs
  chrs <- fread("/scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.daphnids.delim", header=F)
  chrs <- unique(chrs$V1)
  chrs

### fai
  fai <- fread(list.files("/scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/", "fai", full.name=T)[1])
  setkey(fai, V1)

  use <- fai[J(chrs)]

### make windows

  wins <- foreach(i=1:dim(use)[1], .combine="rbind")%do%{
    data.table(chr=use[i]$V1,
                start=seq(from=1, to=use[i]$V2, by=step.bp),
                stop=seq(from=1, to=use[i]$V2, by=step.bp)+window.bp)
  }

### add special windows
  load("/project/berglandlab/alan/gprime_peaks.replicates.250K.05.Rdata")

  wins <- rbind(wins,
                data.table(chr=peaks$CHROM, start=peaks$posMaxGprime-(window.bp/2), stop=peaks$posMaxGprime+(window.bp/2)))


  write.table(wins, file="/scratch/aob2x/daphnia_hwe_sims/popPhase/windows.delim", quote=F, row.names=F, col.names=F, sep=",")
