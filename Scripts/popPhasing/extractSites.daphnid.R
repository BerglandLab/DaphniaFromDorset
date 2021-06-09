#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)

### load
  load("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/dpfiltsnps_20200623.Rdata")

### export
  write.table(dpfiltsnps[,c("chr", "pos"), with=F], file="/scratch/aob2x/daphnia_hwe_sims/popPhase/daphnid.sites", quote=F, row.names=F, col.names=F)
