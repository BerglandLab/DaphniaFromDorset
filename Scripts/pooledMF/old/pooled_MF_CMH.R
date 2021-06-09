### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

## install packges
  #install.packages("stringi")
  #install.packages("matrixStats")
  #install.packages("~/poolSeq.tar.gz", repos=NULL, type="source")

### libraries
  library(poolSeq)
  library(data.table)
  library(SN)
### load data
  #load("/project/berglandlab/Karen/SingleMomsMales20182019/PooledMomsMales/totrdfilt.Rdata")
  load("/nv/vol186/bergland-lab/alan/totalADRDlongall.Rdata")

  geno[,effRD:=floor(RRD)]
  geno[,effPA:=round(propalt*effRD)/effRD]

### CMH test
    ### D8 males
      pond.i<-"D8"


      #D8 <- geno[pond==pond.i][,list(A0=RD, a0=AD), list(Sample, chr, pos)]
      D8 <- geno[pond==pond.i][,list(A0=effRD*(1-effPA), a0=effRD*(effPA)), list(Sample, chr, pos)]

      A.tmp <- as.matrix(dcast(D8[,-"a0", with=F], Sample ~ chr+pos)[,-1])
      a.tmp <- as.matrix(dcast(D8[,-"A0", with=F], Sample ~ chr+pos)[,-1])

      missing <- apply(A.tmp, 2, function(x) any(is.na(x)))

      A.tmp <- A.tmp[,!missing]
      a.tmp <- a.tmp[,!missing]

      A0 <- A.tmp[1:2,]
      At <- A.tmp[3:4,]

      a0 <- a.tmp[1:2,]
      at <- a.tmp[3:4,]


      o <- cmh.test(A0=A0, a0=a0, At=At, at=at, log=T, min.cnt=1)


      dt <- data.table(p=o, snp=dimnames(A.tmp)[[2]])
      D8[,snp:=paste(chr, pos, sep="_")]

      setkey(dt, snp)
      setkey(D8, snp)

      dt <- merge(dt, D8)
      dt.ag <- dt[,list(p=mean(p)), list(snp)]
      dt.ag[,q:=p.adjust(10^(-p), "fdr")]
      dt.ag[q<1]


      tab <- dt[,list(T=sum(p>2, na.rm=T), F=sum(p<=2, na.rm=T)), list(chr)]
      chisq.test(], -1])



      load("/mnt/sammas_storage/bergland-lab/alan/cmh_mf.Rdata")


      load("/nv/vol186/bergland-lab/alan/cmh_mf.Rdata")
