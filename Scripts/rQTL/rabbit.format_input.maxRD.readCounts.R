#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### makes rabbit input
  args = commandArgs(trailingOnly=TRUE)
  chr.i <- as.character(args[1])
  maxcM <- as.numeric(args[2])
  f1s.set <- as.character(args[3])
  #chr.i <- "Scaffold_1863_HRSCAF_2081"; maxcM=10; f1s.set <- "all_AxC"

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

  }
  f1s <- data.table(cloneid=f1s)

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
  #snp.dt.ag <- snp.dt[,.N,chr]
  #write.table(snp.dt.ag, file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase/chrs.csv", quote=F, row.names=T, col.names=F, sep=",")

### make large input file

    ### first, find informative sites

      seqSetFilter(genofile,
                  sample.id=c(sc[SC=="A"][which.max(medrd)]$clone,
                              sc[SC=="C"][which.max(medrd)]$clone),
                  variant.id=snp.dt[J(chr.i)][numAlleles==2]$id)

      dosage <- as.data.table(t(seqGetData(genofile, "$dosage")))
      setnames(dosage, names(dosage), seqGetData(genofile, "sample.id"))
      setnames(dosage, sc[SC=="A"][which.max(medrd)]$clone, "A")
      setnames(dosage, sc[SC=="C"][which.max(medrd)]$clone, "C")

      dosage[,id:=snp.dt[J(chr.i)][numAlleles==2]$id]

      dosage[,use:=F]
      dosage[(A==1 & C==0) | (A==1 & C==2) | (A==0 & C==1) | (A==2 & C==1) | (A==1 & C==1), use:=T]

    ### second, pull out depths

      seqSetFilter(genofile,
                  sample.id=c(sc[SC=="A"][which.max(medrd)]$clone,
                              sc[SC=="C"][which.max(medrd)]$clone,
                              f1s$cloneid),
                  variant.id=dosage[use==T]$id)

      alleleDepths <- seqGetData(genofile, "annotation/format/AD")$data
      refDepth <- altDepth[,seq(from=1, to=dim(altDepth)[2]-1, by=2)]
      altDepth <- altDepth[,seq(from=2, to=dim(altDepth)[2], by=2)]

      f1.ord <- data.table(clone=seqGetData(genofile, "sample.id"))
      f1.ord <- merge(f1.ord, sc[,c("clone", "SC"), with=F])

      parents <- foreach(ind.i=c("A", "C"), .combine="rbind")%do%{
        #ind.i<-"A"
        i <- which(f1.ord[SC==ind.i]$clone==seqGetData(genofile, "sample.id"))

        tmp <- paste(refDepth[i,], altDepth[i,], sep="|")
        tmp <- matrix(c(ind.i, tmp), nrow=1)
        tmp
      }

      offspring <- foreach(ind.i=f1.ord[!SC%in%c("A", "C")]$clone, .combine="rbind")%do%{
        #ind.i<-"A"
        i <- which(f1.ord[clone==ind.i]$clone==seqGetData(genofile, "sample.id"))

        tmp <- paste(refDepth[i,], altDepth[i,], sep="|")
        tmp <- matrix(c(ind.i, tmp), nrow=1)
        tmp
      }
      marker <- matrix(c("marker", seqGetData(genofile, "variant.id")), nrow=1)
      #chr <- matrix(c("chromosome", rep(NA, dim(genomat)[1])), nrow=1)
      #pos <- matrix(c("pos(cM)", rep(NA, dim(genomat)[1])), nrow=1)
      chr <- matrix(c("chromosome", rep(as.numeric(as.factor(chr.i)), dim(marker)[2]-1)), nrow=1)
      pos <- matrix(c("chromosome", seq(from=0, to=maxcM, length.out=dim(marker)[2]-1)), nrow=1)

      header <- do.call("rbind", list(marker, chr, pos))

      out <- do.call("rbind", list(header, parents, offspring))


###

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

writeLines( "Pedigree-Information,DesignPedigree\nGeneration,MemberID,Female=1/Male=2/Hermaphrodite=0,MotherID,FatherID\n0,1,1,0,0\n0,2,2,0,0\n1,3,0,1,2\nPedigree-Information,SampleInfor\nProgenyLine,MemberID,Funnelcode",
             con=ped.fn
           )

f1s[,id:=3]
f1s[,fc:="1-2"]

 write.table(f1s,
             file=ped.fn,
             quote=FALSE,
             row.names=FALSE,
             col.names=FALSE,
             sep=",",
             na="NA",
             append=TRUE)
