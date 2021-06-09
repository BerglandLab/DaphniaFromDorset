#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


### libraries
  library(data.table)
  library(foreach)
  library(ape)
  library(doMC)
  registerDoMC(20)

### input files: genome-wide distribution


  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/popPhase/trees250K/", ".info", full.names=T)
  length(fn)

  regions <- foreach(i=fn)%dopar%{
    message(i)
    #i<-fn[1]
    region <- fread(i)
    region[,N:=V2-(V3+V4+V5+V6)]
    region[,window:=tstrsplit(V1, ";")[[2]]]
    region[,list(N=mean(N), size=mean(V2)), window]
  }
  regions <- rbindlist(regions)

### save
  save(regions, file="~/regions_250K.Rdata")

  #save(cdl.genome, cdl.qtl, cdl.o, cdl.tree, cdl.list, file="~/cdlo.Rdata")
