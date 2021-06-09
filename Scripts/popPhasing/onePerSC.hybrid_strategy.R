
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(SeqArray)

### open GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

  snpFilter <- fread("/scratch/aob2x/daphnia_hwe_sims/popPhase/daphnid.sites")

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
  snp.dt <- merge(snpFilter, snp.dt)



### clone metadata file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  sc <- sc[Nonindependent==0 ]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]

### hard filtering of SC
  sc.ag <- sc[LabGenerated==F & Species=="pulex", list(clone=clone[which.max(medrd)][1]), list(SC.uniq, Species)]
  #sc.ag[,clone:=paste(Species, clone, sep="_")]

### get 12 chromosomes
  ### load in filter file
    load("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/dpfiltsnps_20200623.Rdata")

  ### get 12 chrs
    chrs <- unique(dpfiltsnps$chr)

### make long job file, 1 per chr per clone
  sc.ag.chr <- foreach(i=chrs, .combine="rbind")%do%cbind(data.table(chr=i), sc.ag)

### write to disk
  write.table(sc.ag.chr[,c("chr", "clone"),with=F], file="/scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.hybrid_strategy.delim", quote=F, row.names=F, col.names=F, sep="\t")



### get  obtusa and pulicaria allele frequencies

  consensus <- foreach(sc.i=c("obtusa", "pulicaria"), .combine="cbind")%do%{
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sc[Species==sc.i]$clone, variant.id=snp.dt$id)

    data.table(af=seqAlleleFreq(genofile, ref.allele=1L)) ### alternate allele
  }
  setnames(consensus, c(1,2), c("af.obtusa", "af.pulicaria"))
  consensus <- cbind(consensus, snp.dt)


### make consensus call
  consensus[af.obtusa==0, geno.obtusa:='0|0']
  consensus[af.obtusa==1, geno.obtusa:='1|1']
  consensus[round(af.obtusa,2)==0.5, geno.obtusa:=sample(c('0|0', '1|1'), replace=T, size=sum(round(consensus$af.obtusa, 2)==0.5, na.rm=T))]
  consensus[round(af.obtusa,2)>0.5, geno.obtusa:='1|1']
  consensus[round(af.obtusa,2)<0.5, geno.obtusa:='0|0']
  consensus[is.na(geno.obtusa), geno.obtusa:=".|."]

  consensus[af.pulicaria==0, geno.pulicaria:='0|0']
  consensus[af.pulicaria==1, geno.pulicaria:='1|1']
  consensus[round(af.pulicaria,2)==0.5, geno.pulicaria:=sample(c('0|0', '1|1'), replace=T, size=sum(round(consensus$af.pulicaria, 2)==0.5, na.rm=T))]
  consensus[round(af.pulicaria,2)>0.5, geno.pulicaria:='1|1']
  consensus[round(af.pulicaria,2)<0.5, geno.pulicaria:='0|0']
  consensus[is.na(geno.pulicaria), geno.pulicaria:=".|."]

  table(consensus$geno.obtusa)
  table(consensus$geno.pulicaria)

### format VCF
  vcf <- data.table("#CHROM"=consensus$chr, POS=consensus$pos, ID='.',
                    REF=snp.dt$ref, ALT=snp.dt$alt,
                    QUAL=".", FILTER=".", INFO=".", FORMAT="GT",
                    pulicaria=consensus$geno.pulicaria,
                    obtusa=consensus$geno.obtusa)

### extract header
  system("cat /project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.vcf | head -n 10000 | grep '##' > /scratch/aob2x/daphnia_hwe_sims/popPhase/outgroup.vcf")
  system('tail -n10 /scratch/aob2x/daphnia_hwe_sims/popPhase/outgroup.vcf')

  write.table(vcf, file='/scratch/aob2x/daphnia_hwe_sims/popPhase/outgroup.vcf', append=T, quote=F, row.names=F, sep="\t")





#consensus[!is.na(af.obtusa),obtusa.geno := unlist(sapply(consensus[!is.na(af.obtusa)]$af.obtusa, function(x) c("11","12","22")[which.min(abs(x-c(0,.5,1)))]))]
#consensus[!is.na(af.pulicaria),pulicaria.geno := unlist(sapply(consensus[!is.na(af.pulicaria)]$af.pulicaria, function(x) c("11","12","22")[which.min(abs(x-c(0,.5,1)))]))]

con
