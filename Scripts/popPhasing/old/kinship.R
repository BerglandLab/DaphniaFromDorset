#ijob -c1 -p standard -A berglandlab
#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(SNPRelate)

### load SuperClone & SNP filter file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/CloneInfoFilePulexandObtusa_withmedrd_20200207")
  snps2use <- fread("/project/berglandlab/Karen/MappingDec2019/snpsvarpulexpresentinhalf_table_20200207")
  setnames(snps2use, "variant.ids", "id")
  snps2use[,use:=T]

  setkey(snps2use, chr, pos)

  sc[,SC.uniq:=paste(SC, SCnum, sep="_")]

### open GDS

### open GDS
  #genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"), pos=seqGetData(genofile, "position"), id=seqGetData(genofile, "variant.id"), key="id")
  setkey(snp.dt, chr, pos)

  snp.dt <- merge(snp.dt, snps2use)

  snp.dt[,use.chr:=F]
  snp.dt[chr%in%snp.dt[,.N,chr][N>1000]$chr, use.chr:=T]
