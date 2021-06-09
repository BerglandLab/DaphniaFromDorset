#ijob -c1 -p standard -A berglandlab
#module load intel/18.0 intelmpi/18.0 R/3.6.3; R



### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### open GDS
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

### load in filter file
  snpFilter <- fread("snpsvarpulexpresentinhalf_table_20200623")

### make snp.dt
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       key="chr")
  setkey(snpFilter, chr, pos)
  setkey(snp.dt, chr, pos)

  snp.dt <- merge(snpFilter, snp.dt)

### 1. identify informative differences between A & C
  ac.fd <- foreach(sc.i=c("A", "C"), .combine="cbind")%do%{
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sc[SC==sc.i]$clone, variant.id=snp.dt$id)

    data.table(af=seqAlleleFreq(genofile, ref.allele=1L)) ### alternate allele
  }
  setnames(ac.fd, c(1,2), c("af.A", "af.C"))
  ac.fd <- cbind(ac.fd, snp.dt)


  ac.fd[!is.na(af.A),A.geno := unlist(sapply(ac.fd[!is.na(af.A)]$af.A, function(x) c("11","12","22")[which.min(abs(x-c(0,.5,1)))]))]
  ac.fd[!is.na(af.C),C.geno := unlist(sapply(ac.fd[!is.na(af.C)]$af.C, function(x) c("11","12","22")[which.min(abs(x-c(0,.5,1)))]))]

  ac.fd[!is.na(af.A),A.delta := unlist(sapply(ac.fd[!is.na(af.A)]$af.A, function(x) min(abs(x-c(0,.5,1)))))]
  ac.fd[!is.na(af.C),C.delta := unlist(sapply(ac.fd[!is.na(af.C)]$af.C, function(x) min(abs(x-c(0,.5,1)))))]

  ac.fd <- ac.fd[A.delta < 0.05 & C.delta < 0.05]

  ac.inform <- ac.fd[(A.geno=="12" & C.geno=="11") |
                     (A.geno=="12" & C.geno=="22") |
                     (A.geno=="11" & C.geno=="12") |
                     (A.geno=="22" & C.geno=="12") |
                     (A.geno=="12" & C.geno=="12") ]

### 2. Identify F1s to use
  f1s <- sc[AxCF1Hybrid==1]$clone

### 3. make majority rule F1s, A & B
  setkey(sc, clone)
  sc.f1.ag <- sc[J(f1s)][SC!="OO" & SC!="AxCF1",.N,list(SC)]
  sc.f1.ag <- rbind(sc.f1.ag, data.table(N=100, SC=c("A", "C")))[order(N, decreasing=T)]

  nF1s <- 4

  f1.cons <- foreach(i=1:(2+nF1s), .combine="cbind")%do%{
    #i<-1
    print(i)
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sc[SC==sc.f1.ag$SC[i]]$clone, variant.id=ac.inform$id)

    tmp <- data.table(af=seqAlleleFreq(genofile, ref.allele=1L))
    #tmp <- cbind(tmp, snp.dt)

    tmp[!is.na(af), geno := unlist(sapply(tmp[!is.na(af)]$af, function(x) c("0/0","0/1","1/1")[which.min(abs(x-c(0,.5,1)))]))]

    tmp[,"geno",with=F]
  }
  setnames(f1.cons, c(1:dim(f1.cons)[2]), sc.f1.ag[1:(2+nF1s)]$SC)

  proto.vcf <- cbind(ac.inform, f1.cons)

  vcf <- data.table('#CHROM'=ac.inform$chr,
                    POS=ac.inform$pos,
                    ID=paste("snp", ac.inform$id, sep="_"),
                    REF=seqGetData(genofile, "$ref"),
                    ALT=seqGetData(genofile, "$alt"),
                    QUAL=".",
                    FILTER="PASS",
                    INFO=".",
                    FORMAT="GT")
  vcf <- cbind(vcf, f1.cons)

### 4. write VCF file & PED file
  seqSetFilter(genofile, variant.id=1)
  seqGDS2VCF(genofile,
             "/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.vcf",
             info.var=character(0), fmt.var=character(0),
             verbose=TRUE)

  system("grep '##' /scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.vcf > /scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.vcf")

  write.table(vcf, file="/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.vcf", sep="\t", append=T, col.names=T, quote=F, row.names=F)

  ped <- data.table(fam="family1",
                    iid=names(f1.cons)[-c(1,2)],
                    pid="C",
                    mid="A",
                    foo="N", bar="A")

  write.table(ped, sep="\t", quote=F, row.names=F, col.names=F, file="/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.ped")













  het.count.ag <- het.count[,list(RR=c(mean(fRR[A.geno==0]), mean(fRR[B.geno==0])),
                                  RA=c(mean(fRA[A.geno==1]), mean(fRA[B.geno==1])),
                                  AA=c(mean(fAA[A.geno==2]), mean(fAA[B.geno==2])),
                                  n=sum(n),
                                  SC=c("A", "B")),
                             list(sample.id)]
  het.count.ag[,mu:=(RR+RA+AA)/3]


### make small, dummy vcf
    seqResetFilter(genofile)

    nF1 <- 2
    seqSetFilter(genofile, sample.id=c(het.count.ag[SC=="A"][which.max(mu)]$sample.id,
                                       het.count.ag[SC=="B"][which.max(mu)]$sample.id,
                                       f1.set[1:nF1]),
                           variant.id=snp.dt[use.chr==T]$id)

     seqGDS2VCF(genofile, vcf.fn="/scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.10f1.vcf.gz", info.var=character(0), fmt.var=character(0))

     ped <- data.table(fam="family1",
                       iid=f1.set[1:nF1],
                       pid=het.count.ag[SC=="B"][which.max(mu)]$sample.id,
                       mid=het.count.ag[SC=="A"][which.max(mu)]$sample.id,
                       foo="N", bar="A")

     write.table(ped, sep="\t", quote=F, row.names=F, col.names=F, file="/scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.10f1.ped")

### make ped file
    ped <- data.table(fam="family1",
                      iid=f1.set,
                      pid=sc[SC=="A"]$clone[1],
                      mid=sc[SC=="B"]$clone[1],
                      foo="N", bar="A")

    write.table(ped, sep="\t", quote=F, row.names=F, col.names=F, file="/scratch/aob2x/daphnia_hwe_sims/trioPhase/allTrio.ped")





  save(het.count, file="/nv/vol186/bergland-lab/alan/f1q.Rdata")


  library(data.table)
  library(ggplot2)

  load("/mnt/sammas_storage/bergland-lab/alan/f1q.Rdata")

  hist(het.count[A.geno==0 & B.geno==2]$fRA)
  abline(v=.90, col="red")


  het.count[sample.id%in%het.count[n>10][A.geno==0 & B.geno==2 & fRA>.9]$sample.id]
