#ijob -c1 -p standard -A berglandlab
#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone=
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### load SNP filter file
  load("snpsvarpulexpresentinhalf_20200623.Rdata")

### open GDS
  genofile <- seqOpen("MapJune2020_ann.seq.gds")
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"), pos=seqGetData(genofile, "position"), id=seqGetData(genofile, "variant.id"), key="id")
  use <- data.table(id=snpsvarpulexpresentinhalf, pass=T, key="id")
  snp.dt <- merge(snp.dt, use)

  snp.dt[,use.chr:=F]
  snp.dt[chr%in%snp.dt[,.N,chr][N>1000]$chr, use.chr:=T]

### 1. identify fixed difference between A&C
  ac.fd <- foreach(sc.i=c("A", "C"), .combine="cbind")%do%{
    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sc[SC==sc.i]$clone, variant.id=snp.dt$id)

    data.table(af=seqAlleleFreq(genofile))
  }
  setnames(ac.fd, c(1,2), c("af.A", "af.C"))
  ac.fd <- cbind(ac.fd, snp.dt)
  ac.fd[!is.na(af.A),A.geno := unlist(sapply(ac.fd[!is.na(af.A)]$af.A, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  ac.fd[!is.na(af.C),C.geno := unlist(sapply(ac.fd[!is.na(af.C)]$af.C, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]

### 2. Identify informative sites
  ac.inform <- ac.fd[(A.geno==0 & C.geno==2) | (A.geno==1 & C.geno==0) | (A.geno==0 & C.geno==1) | (A.geno==2 & C.geno==0)]


  ### 3. make majority rule F1s, A & B
    setkey(sc, clone)
    sc.f1.ag <- sc[J(f1.set)][SC!="OO" & SC!="AxCF1",.N,list(SC)]
    sc.f1.ag <- rbind(sc.f1.ag, data.table(N=100, SC=c("A", "C")))[order(N, decreasing=T)]

    nF1s <- 4

    f1.cons <- foreach(i=1:(2+nF1s), .combine="cbind")%do%{
      #i<-1
      print(i)
      seqResetFilter(genofile)
      seqSetFilter(genofile, sample.id=sc[SC==sc.f1.ag$SC[i]]$clone, variant.id=ac.inform$id)

      tmp <- data.table(af=seqAlleleFreq(genofile))
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
    seqGDS2VCF(genofile, "/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.vcf", info.var=character(0), fmt.var=character(0),
      verbose=TRUE)
    system("grep '##' /scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.vcf > /scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.vcf")

    write.table(vcf, file="/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.vcf", sep="\t", append=T, col.names=T, quote=F, row.names=F)

    ped <- data.table(fam="family1",
                      iid=names(f1.cons)[-c(1,2)],
                      pid="C",
                      mid="A",
                      foo="N", bar="A")

    write.table(ped, sep="\t", quote=F, row.names=F, col.names=F, file="/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.ped")
