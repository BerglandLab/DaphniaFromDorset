#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### makes rabbit input
  args = commandArgs(trailingOnly=TRUE)
  chr.i <- as.character(args[1])
  maxcM <- as.numeric(args[2])
  f1s.set <- as.character(args[3])
  #chr.i <- "Scaffold_1863_HRSCAF_2081"; maxcM=10; f1s.set <- "all"

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### which F1s?
  #f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.onlyPheno.delim")
  #f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.allF1s.delim")
  #f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.all_AxC_F1s.delim")


  if(f1s.set=="onlyPheno_AxC") {
    f1s <- sc[AxCF1Hybrid==1][OneLiterPheno==1]$clone

  } else if (f1s.set=="wildF1s_AxC"){
    f1s <- sc[AxCF1Hybrid==1][OneLiterPheno==0]$clone

  } else if(f1s.set=="all_AxC") {
    f1s <- sc[AxCF1Hybrid==1]$clone

  } else if(f1s.set=="all_CxC") {
    f1s <- sc[OneLiterPheno==1][AxCF1Hybrid==0][SC=="selfedC"]$clone

  } else if(f1s.set=="all") {
    f1s <- c(sc[AxCF1Hybrid==1]$clone,
             sc[OneLiterPheno==1][AxCF1Hybrid==0][SC=="selfedC"]$clone)
  }
  f1s <- data.table(clone=f1s)


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



### make majority rule (consensus) genotype calls for
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
                     (A.geno=="12" & C.geno=="12") |
                     (A.geno=="11" & C.geno=="22") |
                     (A.geno=="22" & C.geno=="11") ]

  ac.inform <- ac.inform[chr==chr.i]
  ac.inform[,pos.bin:=round(pos/1e5)]


### select sites in F1s with lowest amount of missing data
  seqResetFilter(genofile)
  seqSetFilter(genofile,
              sample.id=f1s$clone,
              variant.id=ac.inform$id)

  mr <- data.table(id=seqGetData(genofile, "variant.id"),
                    mr=seqMissing(genofile))
  mr <- merge(mr, snp.dt, by="id")


  ### check to see if missing rates are homogeneously distributed throughout the genome.
    mr[,pos.bin:=round(pos/1e5)]
    mr.ag <- mr[,list(nLow=sum(mr<.25), n=length(mr)), list(pos.bin)]
    summary(mr.ag$nLow/mr.ag$n)

  ### trim out position bins with high rates of missing data (i.e., when nLow/n is high. mr=missing rate so a low is good; we want windows with a lot of sites with low missing rates, i.e. with y>.5)

    #ggplot(data=mr.ag[nLow/n>.5], aes(y=nLow/n, x=pos.bin)) + geom_line()

  ### select windows with low rates of high missing rat ( 50% )
      setkey(ac.inform, pos.bin)
      ac.inform <- ac.inform[J(mr.ag[nLow/n >.5]$pos.bin)]

  ### trim to sites with low rates of missing data
      ac.inform <- merge(ac.inform, mr[,c("id", "mr"), with=F], by="id")
      ac.inform <- ac.inform[mr<.25]


  ### subsample
    ac.inform.ag <- ac.inform[,list(n=length(id)), list(pos.bin)]

    sample.fun <- function(x, n) {
      if(length(x)<=n) return(x)
      if(length(x)>n) return(sort(as.integer(sample(as.character(x), size=n))))
    }

    set.seed(1234)

    ac.inform.sub <- ac.inform[,list(id=sample.fun(id, 50000)), list(pos.bin)]

    setkey(ac.inform.sub, pos.bin, id)
    setkey(ac.inform, pos.bin, id)
    ac.inform <- merge(ac.inform, ac.inform.sub)


    ac.inform[,list(n=length(id)), list(pos.bin)]

  ### make parents
   A.parent <- c("A", ac.inform$A.geno)
   C.parent <- c("C", ac.inform$C.geno)

   head(A.parent)
   head(C.parent)
   parents <- rbind(A.parent, C.parent)
   head(parents[,1:10])

### load & format offspring
  seqResetFilter(genofile)
  seqSetFilter(genofile,
              sample.id=f1s$clone,
              variant.id=ac.inform$id)


  genomat <- as.data.table(t(seqGetData(genofile, "$dosage")))
  setnames(genomat, seqGetData(genofile, "sample.id"))

  genomat[,id:=seqGetData(genofile, "variant.id")]

### check
  genomat.l <- melt(genomat, id.vars="id")
  genomat.l.ag <- genomat.l[,list(n22=sum(value==0, na.rm=T), n12=sum(value==1, na.rm=T), n11=sum(value==2, na.rm=T)), list(id)]

  gp <- merge(genomat.l.ag, ac.inform, by="id")





  offspring <- foreach(ind.i=f1s$clone, .combine="rbind", .errorhandling="remove")%do%{
    tmp <- t(as.matrix(genomat[,ind.i, with=F]))
    tmp[tmp=="0"] <- "2N"
    #tmp[tmp=="1"] <- sample(c("1N","2N"), dim(tmp)[1], replace=T)
    tmp[tmp=="1"] <- "12"
    tmp[tmp=="2"] <- "1N"
    tmp[is.na(tmp)] <- "NN"
    cbind(matrix(ind.i, ncol=1), tmp)
  }


  dim(offspring)
  offspring[1:5,1:10]

### make header

  marker <- matrix(c("marker", seqGetData(genofile, "variant.id")), nrow=1)
  #chr <- matrix(c("chromosome", rep(NA, dim(genomat)[1])), nrow=1)
  #pos <- matrix(c("pos(cM)", rep(NA, dim(genomat)[1])), nrow=1)
  chr <- matrix(c("chromosome", rep(as.numeric(as.factor(chr.i)), dim(marker)[2]-1)), nrow=1)
  pos <- matrix(c("pos(cM)", seq(from=0, to=maxcM, length.out=dim(marker)[2]-1)), nrow=1)

  header <- do.call("rbind", list(marker, chr, pos))


### combine
  out <- do.call("rbind", list(header, parents, offspring))
  rownames(out) <- NULL
  out[1:7,1:4]


### write

  out.fn <- paste("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_", maxcM, "cm/", chr.i, "/", chr.i, ".all.in", sep="")

  writeLines( paste("#founders,",2, sep=""),
               con=out.fn
             )
  options(scipen=999)

   write.table(out,
               file=out.fn,
               quote=FALSE,
               row.names=FALSE,
               col.names=FALSE,
               sep=",",
               na="NA",
               append=TRUE)

### make ped file
  ped.fn <- paste("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_", maxcM, "cm/", chr.i, "/", chr.i, ".ped", sep="")

  if(f1s.set!="all") {
    writeLines( "Pedigree-Information,DesignPedigree\nGeneration,MemberID,Female=1/Male=2/Hermaphrodite=0,MotherID,FatherID\n0,1,1,0,0\n0,2,2,0,0\n1,3,0,1,2\nPedigree-Information,SampleInfor\nProgenyLine,MemberID,Funnelcode",
                 con=ped.fn
               )

    f1s[,id:=3]
    f1s[,fc:="1-2"]
  } else if(f1s.set=="all" ) {
    writeLines( "Pedigree-Information,DesignPedigree\nGeneration,MemberID,Female=1/Male=2/Hermaphrodite=0,MotherID,FatherID\n0,1,1,0,0\n0,2,0,0,0\n1,3,0,1,2\n1,4,0,2,2\nPedigree-Information,SampleInfor\nProgenyLine,MemberID,Funnelcode",
                 con=ped.fn
               )

    f1s[clone%in%sc[AxCF1Hybrid==1]$clone, id:=3]
    f1s[clone%in%sc[OneLiterPheno==1][AxCF1Hybrid==0][SC=="selfedC"]$clone, id:=4]
    f1s[,fc:="1-2"]


  }
   write.table(f1s,
               file=ped.fn,
               quote=FALSE,
               row.names=FALSE,
               col.names=FALSE,
               sep=",",
               na="NA",
               append=TRUE)
