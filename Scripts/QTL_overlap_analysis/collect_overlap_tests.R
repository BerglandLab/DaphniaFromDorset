# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

#permUse <-0; crossType<-"CxC"
### libraries

  library(data.table)
  library(foreach)

### load
  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/overlap_test", "overlap", full.name=T)

  o <- foreach(fn.i=fn)%do%{
    #fn.i <- fn[2]

    tmp <- fread(fn.i)
    tmp[,cross:=ifelse(grepl("AxC", fn.i), "AxC", "CxC")]
    tmp
  }
  overlap.perm <- rbindlist(o, fill=T)

### save
  save(overlap.perm, file="~/overlap_perm.Rdata")


### scp aob2x@rivanna.hpc.virginia.edu:~/overlap_perm.Rdata ~/.
  library(data.table)
  overlap.perm[,list(pr=mean(z[perm==0] > z[perm!=0], na.rm=T)), list(pheno, cross)]
