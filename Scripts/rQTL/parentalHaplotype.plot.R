# scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/parental.Rdata ~/.


### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(patchwork)

### load data
  load("~/peaks.Rdata")
  load("~/parental.Rdata")
  setkey(parental, chr, pos)

### function
  plotParentalHaplos <- function(chr.i, maxPos, window=10000) {
    #chr.i=peaks[which.max(maxGprime)]$CHROM; start=peaks[which.max(maxGprime)]$posMaxGprime - 5000; stop=peaks[which.max(maxGprime)]$posMaxGprime + 5000

    start<-maxPos-window
    stop<-maxPos+window
    tmp <- parental[J(data.table(chr=chr.i, pos=start:stop)), nomatch=0]
    tmp[,id.x:=rank(id, ties="min")]

    ggplot(data=tmp, aes(x=id.x, y=allele, fill=as.factor(value))) + geom_tile()

  }

  i<-12
  plotParentalHaplos(chr.i=peaks$CHROM[i],
                    maxPos=peaks$posMaxGprime[i],
                    window=10000)
