#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(SeqArray)
  options(scipen=999)
### open GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

  snpFilter <- fread("/scratch/aob2x/daphnia_hwe_sims/popPhase/daphnid.sites")
  snpFilter[,use:=T]
  setnames(snpFilter, c("V1", "V2"), c("chr", "pos"))
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       key="chr",
                       ref=seqGetData(genofile, "$ref"),
                       alt=seqGetData(genofile, "$alt"))

  setkey(snpFilter, chr, pos)
  setkey(snp.dt, chr, pos)
  snp.dt
  snp.dt <- merge(snpFilter, snp.dt, all.y=T)

  snp.dt[is.na(use)]

  bed <- snp.dt[is.na(use)][,list(chr=chr, start=pos-1, stop=pos)]
  setkey(bed, chr, start)

  write.table(bed, file="/scratch/aob2x/daphnia_hwe_sims/popPhase/badSites.bed", quote=F, row.names=F, col.names=F)
  
