#module load intel/18.0 intelmpi/18.0 R/3.6.0; R

### libraries
  library(data.table)
  #library(GWASTools)
  library(SeqArray)
  library(SNPRelate)
  source("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/ibdAssignRelatedness.R")
  #library(sequoia)

### open GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")

### clone metadata file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/Superclones201617182019withObtusaandPulicaria_kingcorr_20200402_wmedrd.txt")

  sc <- sc[Nonindependent==0]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]

### tack in missing data
  sc.missing <- data.table(clone=seqGetData(genofile, "sample.id"), missing=seqMissing(genofile, per.variant=F))
  sc <- merge(sc, sc.missing, by="clone")

### hard filtering of SC
  #sc <- sc[SC!="OO"]

  # clone[which.max(medrd)][1]
  sc.ag <- sc[Species=="pulex" & !grepl("W|D10", clone) & grepl("D8|DBunk", clone) & missing<.15,
              list(clone=clone[which.min(missing)][1],
                   year=year[which.min(year)][1]),
              list(SC.uniq)]

  sc.ag[,age:=max(sc.ag$year) - year + 1]
  table(sc.ag$SC.uniq)
  table(tstrsplit(sc.ag$clone, "_")[[3]])


### filtered SNP set, LD-pruned by Karen
  snps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/finalsetsnpset01pulex_table_wpulicaria_20200401")
  snps.ag <- snps[,list(.N), chr]

  snps[,goodChr:=F]
  snps[chr%in%snps.ag[N>5000]$chr, goodChr:=T]

### set filter
  seqResetFilter(genofile)
  seqSetFilterPos(genofile, chr=snps[goodChr==T]$chr, pos=snps[goodChr==T]$pos)

  snps.use <- data.table(chr=seqGetData(genofile, "chromosome"),
                           pos=seqGetData(genofile, "position"),
                          id=seqGetData(genofile, "variant.id"))

  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=snps.use$id, sample.id=sc.ag$clone)

  snps.use <- data.table(chr=seqGetData(genofile, "chromosome"),
                           pos=seqGetData(genofile, "position"),
                          id=seqGetData(genofile, "variant.id"),
                          af=seqAlleleFreq(genofile),
                          missing=seqMissing(genofile, per.variant=T))

### export plink files
#seqResetFilter(genofile)
#seqSetFilter(genofile, variant.id=snps.use$id, sample.id=sc.ag$clone)

#seqGDS2VCF(genofile,
#          vcf.fn="/scratch/aob2x/daphnia_hwe_sims/pedigree/onePerSC.LDprune.bcf",
#          info.var=character(), fmt.var=character(), use_Rsamtools=T, verbose=TRUE)


### define inference SNPs
  snps.use[,use:=F]
  snps.use[af>.3 & af<.7 & missing<.05, use:=T]


### Identify "true" F1s
  seqSetFilter(genofile, sample.id=sc.ag[SC.uniq%in%c("A", "C")]$clone, variant.id=snps.use[use==T]$id)
  ab.genomat <- seqGetData(genofile, "$dosage")
  dimnames(ab.genomat) <- list(c("A", "C"), paste("v", snps.use[use==T]$id, sep=""))
  ab.dt <- as.data.table(t(ab.genomat))
  ab.dt[,variant.id:=snps.use[use==T]$id]

  ab.fixed <- #ab.dt[(A==0 & C==2) | (A==2 & C==0)]
  ab.fixed <- ab.dt

### 2. Identify F1s
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=ab.fixed$variant.id, sample.id=sc.ag$clone)

  genomat <- seqGetData(genofile, "$dosage")

  het.count <- foreach(i=1:dim(genomat)[1], .combine="rbind")%do%{
    print(paste(i, dim(genomat)[1], sep=" / "))

    tmp <- ab.fixed
    tmp[,geno:=genomat[i,]]
    nLoci <- dim(tmp[!is.na(A) & !is.na(C) & !is.na(geno)])[1]

    tmp.ag <- tmp[!is.na(A) & !is.na(C),list(nRR=sum(geno==2, na.rm=T),
                                             nRA=sum(geno==1, na.rm=T),
                                             nAA=sum(geno==0, na.rm=T),
                                               n=sum(!is.na(geno))),
                                        list(A.geno=abs(A-2), C.geno=abs(C-2))]
    tmp.ag[,fRR:=nRR/n]
    tmp.ag[,fRA:=nRA/n]
    tmp.ag[,fAA:=nAA/n]

    tmp.ag[,sample.id:=seqGetData(genofile, "sample.id")[i]]
    tmp.ag
  }

  f1.set <- het.count[A.geno==0 & C.geno==2 & fRA>.9]$sample.id
  parent.set <- sc.ag[SC.uniq%in%c("A", "C")]$clone

  family.set <- c(f1.set, parent.set)


### IBD MoM
  ibd.mom <- snpgdsIBDMoM(genofile,
                snp.id=snps.use[use==T]$id,
                sample.id=family.set,
                autosome.only=FALSE,
                maf=.01, verbose=T)


  mom.dt <- foreach(d=grep("k", names(ibd.jac)), .combine="rbind")%do%{
    #d <- 1
    tmp <- ibd.mom[[d]]
    tmp[upper.tri(tmp)] <- NA

    dimnames(tmp) <- list(ibd.jac$sample.id, ibd.jac$sample.id)

    tmp.dt <- as.data.table(melt(tmp))
    setnames(tmp.dt, c("Var1", "Var2", "value"), c("ID1", "ID2", "est"))
    tmp.dt[,coef:=names(ibd.jac)[d]]

    tmp.dt <- na.omit(tmp.dt[ID1!=ID2])

    tmp.dt

  }

  jac.dt[,list(mu=median(est)), list(coef)]

### King
  ### get King Relatedness estimate
    king <- snpgdsIBDKING(genofile,
                  snp.id=snps.use[use==T]$id,
                  sample.id=family.set,
                  autosome.only=FALSE,
                  type=c("KING-robust"))

    #king.cut <- snpgdsIBDSelection(king, kinship.cutoff=1/32)
    king.cut <- king

  ### relatedness
    rel.king <- ibdAssignRelatednessKing(ibs0=king.cut$IBS0, kc=king.cut$kinship, cut.kc.dup=1/(2^(3/2)),
                           cut.kc.fs=1/(2^(5/2)), cut.kc.deg2=1/(2^(7/2)),
                           cut.kc.deg3=1/(2^(9/2)), cut.ibs0.err=0.01)
    table(rel.king)

  ### get IBS values
    ibs <- snpgdsIBS(genofile,
                  snp.id=snps.use[af>.3 & af<.7][missing<.05]$id,
                  sample.id=sc.ag$clone,
                  autosome.only=FALSE)

    ibs.ibs <- ibs$ibs
    ibs.ibs[upper.tri(ibs.ibs)] <- NA
    dimnames(ibs.ibs) <- list(ibs$sample.id, ibs$sample.id)

    ibs.dt <- as.data.table(melt(ibs.ibs))
    setnames(ibs.dt, c("Var1", "Var2", "value"), c("ID1", "ID2", "ibs"))

    ibs.dt[,ID1:=as.character(ID1)]
    ibs.dt[,ID2:=as.character(ID2)]

  ### format IBD & kinship table
    # king
      king.rel <- king$kinship
      king.rel[upper.tri(king.rel)] <- NA
      dimnames(king.rel) <- list(king$sample.id, king$sample.id)

      king.dt <- as.data.table(melt(king.rel))
      setnames(king.dt, c("Var1", "Var2", "value"), c("ID1", "ID2", "kinship"))

      king.dt[,rel:=rel.king]

      king.dt <- king.dt[ID1!=ID2][!is.na(kinship)]
      king.dt[,i:=1:dim(king.dt)[1]]
      king.dt[,ID1:=as.character(ID1)]
      king.dt[,ID2:=as.character(ID2)]

      king.dt <- merge(king.dt, sc.ag, by.x="ID1", by.y="clone")
      king.dt <- merge(king.dt, sc.ag, by.x="ID2", by.y="clone")


      setkey(king.dt, ID1, ID2)
      setkey(ibs.dt, ID1, ID2)

      king.dt <- merge(king.dt, ibs.dt)

      ### define "true" relationship
        king.dt[,true_rel:="FS"]
        king.dt[SC.uniq.x=="A" & SC.uniq.y=="C", true_rel:="U"]
        king.dt[SC.uniq.x=="A" & SC.uniq.y!="C", true_rel:="PO"]
        king.dt[SC.uniq.y=="A", true_rel:="PO"]
        king.dt[SC.uniq.y=="C" & SC.uniq.y!="A", true_rel:="PO"]
        king.dt[SC.uniq.y=="C" & SC.uniq.x=="A", true_rel:="U"]

    ### fold in Jacquard's D7
      jac.dt[coef=="D7"]

      setkey(king.dt, ID1, ID2)
      setkey(jac.dt, ID1, ID2)
      m <- merge(king.dt, jac.dt[coef=="D7"])







      ### ground truth with A x C offspring (C is new B)


        f1s <- unique(king.dt[ID1%in%intersect(king.dt[rel=="PO"][SC.uniq.y=="C"]$ID1,
                                 king.dt[rel=="PO"][SC.uniq.y=="A"]$ID1)][rel=="PO"][SC.uniq.y%in%c("A", "C")]$ID1)

      ### are all of these individuals full sibs?
        foreach(i=1:(length(f1s)-1)%do%{
          foreach(j=(i+1):length(f1s))%do%{
            #i<-1; j<-2
            tmp1 <- f1s[i]
            tmp2 <- f1s[j]

            king.dt[ID1==tmp1][ID2==tmp2]
d
             & ID1==f1s[j]]

          }
        }


      parents <- unique(king.dt[SC.uniq.y%in%c("A", "C")]$ID2)

      king.dt[rel=="Dup"]


  ### Sequoia analysis
    ### extract genomatrix
      seqSetFilter(genofile, sample.id=c(f1s, parents), variant.id=snps.use[af>.3 & af<.7][missing<.01]$id)
      genomat <- seqGetData(genofile, "$dosage")


      genomat <- as.matrix(genomat)
      rownames(genomat) <- gsub("_", "-", seqGetData(genofile, "sample.id"))

      dim(genomat)
      genomat[is.na(genomat)] <- -9
      CheckGeno(genomat)

    ### format age data

      LhD <- data.frame(ID=sc.ag[clone%in%c(f1s, parents)]$clone, Sex=4, BirthYear=sc.ag[clone%in%c(f1s, parents)]$year)
      LhD$ID <- gsub("_", "-", LhD$ID)

    ### save
      save(genomat, LhD, sc.ag, file="/nv/vol186/bergland-lab/alan/seqIn.Rdata")





### IBD
  ### get IBD values
    ibd.mom <- snpgdsIBDMoM(genofile,
                  snp.id=snps.use$id,
                  sample.id=sc.ag$clone,
                  autosome.only=FALSE,
                  maf=.01)

    ibd.jac <- snpgdsIBDMLE(genofile,
                  snp.id=snps.use$id,
                  sample.id=sc.ag$clone,
                  autosome.only=FALSE,
                  maf=.01, method="Jacquard", verbose=T)


  ### save
    save(king, ibd.mom, ibd.jac, file="/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/snprelate_ibd.Rdata")
    load(file="/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/snprelate_ibd.Rdata")
    save(king, ibd.mom, ibd.jac, rel.king, file="/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/snprelate_ibd.Rdata")


### Sequoia analysis
  ### extract genomatrix
    seqSetFilter(genofile, sample.id=sc.ag$clone, variant.id=snps.use[af>.2 & af<.8][missing<.01]$id)
    genomat <- seqGetData(genofile, "$dosage")


    genomat <- as.matrix(genomat)
    rownames(genomat) <- gsub("_", "-", seqGetData(genofile, "sample.id"))

    dim(genomat)
    genomat[is.na(genomat)] <- -9
    CheckGeno(genomat)

  ### format age data
    sc.ag[,clone:=gsub("_", "-", clone)]

    LhD <- data.frame(ID=sc.ag$clone, Sex=4, BirthYear=sc.ag$year)

  ### save
    save(genomat, LhD, file="/nv/vol186/bergland-lab/alan/seqIn.Rdata")

  ###
    library(sequoia)
    load("/mnt/sammas_storage/bergland-lab/alan/seqIn.Rdata")
    LhD$Sex <- 3
    LhD$Sex[LhD$BirthYear!="2018"] <- c(1,2)
    #LhD$BirthYear <- -1

    test <- c("March20-2018-D8-11", "March20-2018-D8-12", "April-2017-D8-151", "April-2017-D8-213")
    test.geno <- genomat[dimnames(genomat)[[1]]%in%test,]
    test.LhD <- LhD[LhD$ID%in%test,]

    test.LhD$Sex <- c(3,3,3,3)

    parOut <- sequoia(GenoM=genomat,
                      LifeHistData=LhD, FindMaybeRel=F,
                      Tfilter = -3, Tassign=.01, MaxMismatch=5, UseAge="yes", CalcLLR=F)

    GetMaybeRel(GenoM=test.geno,
                      LifeHistData=test.LhD

                      , Tfilter = -3, Tassign=.25, MaxMismatch=50, FindMaybeRel=T, UseAge="yes", CalcLLR=T)



                      ,
                      quiet=F, Tfilter = -3, Tassign=.25, MaxMismatch=50, FindMaybeRel=T, UseAge="no", CalcLLR=F)


    test <- c("March20-2018-D8-11", "April-2017-D8-151", "April-2017-D8-213")

    g <- data.table(o.g = genomat[dimnames(genomat)[[1]]==o,],
              p1.g = genomat[dimnames(genomat)[[1]]==p1,],
              p2.g = genomat[dimnames(genomat)[[1]]==p2,])

    table((g$o.g==0 & g$p1.g==2) + (g$o.g==2 & g$p1.g==0))
    table((g$o.g==0 & g$p2.g==2) + (g$o.g==2 & g$p2.g==0))

    8000*1e-04





    ped <- as.data.table(parOut$Pedigree)
    ped[id=="March20-2018-D8-11"]


    mt <- as.data.table(parOut$MaybeTrio)
    mt[id=="March20-2018-D8-11"]
]
    March15-2019-D8-MomPE13   March20-2018-D8-9




    parOut2 <- sequoia(GenoM=genomat, MaxSibIter=10, LifeHistData=LhD, quiet=F, MaxSibshipSize=200, FindMaybeRel=T, MaxMismatch=50)
    save(parOut, file="~/parOut.Rdata")

    maybeRel <- GetMaybeRel(GenoM=genomat)

    ped <- parOut$Pedigree


## cp /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/snprelate_ibd.Rdata /nv/vol186/bergland-lab/alan/snprelate_ibd.Rdata









##### IBD analysis
### libraries
  library(data.table)
  library(FamAgg)

### load
  load("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/snprelate_ibd.Rdata")

### format IBD & kinship table for PRIMUS
  # king
    king.rel <- king$kinship
    king.rel[lower.tri(king.rel)] <- NA
    dimnames(king.rel) <- list(king$sample.id, king$sample.id)

    king.dt <- as.data.table(melt(king.rel))
    setnames(king.dt, c("Var1", "Var2", "value"), c("ID1", "ID2", "kinship"))

    king.dt[,rel:=rel.king]

    king.dt <- king.dt[ID1!=ID2][!is.na(kinship)]
    king.dt[,i:=1:dim(king.dt)[1]]
    king.dt[,ID1:=as.character(ID1)]
    king.dt[,ID2:=as.character(ID2)]


### load truffle data
  truf.ibd <- fread("/scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.vcf.ibd.ibd")


### merge
  setkey(king.dt, ID1, ID2)
  setkey(truf.ibd, ID1, ID2)

  m <- merge(king.dt, truf.ibd)



  save(m, file="~/king_ibd.Rdata")



library(ggplot2)
library(data.table)

load("king_ibd.Rdata")

ggplot(data=m, aes(x=IBD0, y=kinship, color=rel)) + geom_point()






### find 1-degree families
  clones <- unique(c(king.dt$ID1, king.dt$ID2))
  clones.dt <- data.table(individual=clones, father=NA, mother=NA, sibs=NA)
  clones.dt[,individual.n:=as.numeric(as.factor(clones))]

  o <- foreach(c.i=clones.dt$individual, .combine="rbind")%do%{
    #c.i <- clones.dt$individual[1]
    year.i <-



    tmp <- king.dt[ID1==c.i | ID2==c.i]

    PO <- tmp[rel=="PO"]
    FS <- tmp[rel=="FS"]

    data.table(individual=c.i, father=NA, mother=NA,
                            PO=dim(PO)[1],
                            sibs=paste(sort(FS$i), collapse=","),
                            nSibs=length(FS$i))


  }







  # IBD
    # IBD0
      IBD0 <-  ibd.mom$k0
      IBD0[upper.tri(IBD0)] <- NA
      dimnames(IBD0) <- list(ibd.mom$sample.id, ibd.mom$sample.id)

      IBD0.dt <- as.data.table(melt(IBD0))
      setnames(IBD0.dt, c("Var1", "Var2", "value"), c("ID1", "ID2", "ibd0"))

      IBD0.dt <- IBD0.dt[ID1!=ID2][!is.na(ibd0)]

    # IBD1
      IBD1 <-  ibd.mom$k1
      IBD1[upper.tri(IBD1)] <- NA
      dimnames(IBD1) <- list(ibd.mom$sample.id, ibd.mom$sample.id)

      IBD1.dt <- as.data.table(melt(IBD0))
      setnames(IBD1.dt, c("Var1", "Var2", "value"), c("ID1", "ID2", "ibd1"))

      IBD1.dt <- IBD1.dt[ID1!=ID2][!is.na(ibd1)]


### merge
  setkey(king.dt, ID1, ID2)
  setkey(IBD0.dt, ID1, ID2)
  setkey(IBD1.dt, ID1, ID2)

  m <- merge(king.dt, IBD0.dt)
  m <- merge(m, IBD1.dt)
  m[,ibd2:=1 - (m$ibd0 + m$ibd1)]

  m[,pi_hat:= ibd2 + 0.5*ibd1]

### format
  " the input file must have columns 1-4 and 7-10 described in option --plink_ibd, except that the PI_HAT/RELATEDNESS column can be a different measure of relatednes"
  "FID1(1) IID1(2) FID2(3) IID2(4) RT(5) EZ(6) IDB0(7) IBD1(8) IBD2(9) PI_HAT(10)"

  out <- data.table(FID1=m$ID1, IID1=m$ID1, FID2=m$ID2, IID2=m$ID2, RT=0, EZ=0, IBD0=m$ibd0, IBD1=m$ibd1, IBD2=m$ibd2, PI_HAT=m$pi_hat)

  write.table(out, file="/scratch/aob2x/daphnia_hwe_sims/pedigree/snprelatre_output.kinship.delim",
              quote=F, row.names=F, sep="\t")


### write age file
  ages <- data.table(FID=sc.ag$clone, IID=sc.ag$clone, age=sc.ag$age)

  write.table(ages,
            file="/scratch/aob2x/daphnia_hwe_sims/pedigree/snprelatre_output.ages.delim",
                        quote=F, row.names=F, sep="\t")
