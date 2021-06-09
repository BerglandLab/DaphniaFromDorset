### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R



### library
  library(data.table)
  library(foreach)
  library(stringr)

  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out", full.names=T)

  coverage <- foreach(fn.i=fn[1:7], .errorhandling="remove")%do%{
    #fn.i <- fn[1]
    tmp <- fread(paste(fn.i, "QC.summary.txt", sep="/"), header=T)


    tmpp <- tmp[grepl("ReadPairs", FIELD), c("FIELD", "COUNT"), with=F]

    tw <- dcast(tmpp, .~FIELD)

    o <- tw[,list(nogene=ReadPairs_NoGene/(ReadPairs_NoGene + ReadPairs_UniqueGene + ReadPairs_NoGene_Intron + ReadPairs_AmbigGene),
             middle_of_nowhere_rate=ReadPairs_NoGene_MiddleOfNowhere/(ReadPairs_NoGene + ReadPairs_UniqueGene + ReadPairs_NoGene_Intron + ReadPairs_AmbigGene),
             intron_rate=ReadPairs_NoGene_Intron/(ReadPairs_NoGene + ReadPairs_UniqueGene + ReadPairs_NoGene_Intron + ReadPairs_AmbigGene))]
    o[,samp:=last(tstrsplit(fn.i, "/"))]
    melt(o, "samp")


  }
  coverage <- rbindlist(coverage)
  coverage[,list(mu=mean(value)), list(variable)]
  
