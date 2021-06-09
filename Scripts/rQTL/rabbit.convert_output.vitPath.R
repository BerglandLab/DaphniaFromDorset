#module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(SeqArray)
  library(doMC)
  registerDoMC(20)

### make SNP table
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile))
  setkey(snp.dt, id)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### load vit paths
  setwd("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm")

  loadDat <- function(fn) {
    print(fn)
    #fn <- fns[3]

    ### most likely genotype
      pp <- fread(fn)
      pp[,chr:=tstrsplit(V1, " ")[[2]]]
      pp[,clone:=tstrsplit(V1, " ")[[1]]]

    ### return
      return(pp)
  }
  fns <- system("ls /scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/*/*haps*", intern=T)

  ppl <- foreach(x=fns, .errorhandling="remove")%do%loadDat(x)
  ppl <- rbindlist(ppl)


### load parental genotypes

  loadDat <- function(fn) {
    print(fn)
    #fn <- fns[1]
    ### phased & imputed
      pp <- fread(fn, skip=1, nrows=6, header=T, fill=T)
      ppl <- melt(pp, id.vars="marker")

      setnames(ppl, c("marker", "variable"), c("allele", "id"))

      setkey(ppl, "allele")

      ppl <- ppl[J(c("A_Maternal", "A_Paternal", "C_Maternal", "C_Paternal"))]
      ppl[,id:=as.numeric(as.character(id))]

      setkey(ppl, "id")
      setkey(snp.dt, "id")
      ppl <- merge(ppl, snp.dt)

    ### consensus impute

      fni <- gsub(".out.post.csv", ".in",  fn)
      con <- as.matrix(fread(fni, skip=4, nrows=2, header=F, fill=T))
      #con[,1:5]
      id <- as.matrix(fread(fni, skip=1, nrows=1, header=F, fill=T))
      #id[,1:5]
      con.dt <- data.table(id=as.numeric(id[-1]),
                            A=con[1,-1],
                            C=con[2,-1])
      setkey(con.dt, id)
      setkey(ppl, id)
      pplc <- merge(ppl, con.dt)




    return(pplc)
  }

  fns <- system("ls /scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/*/*out.post.csv", intern=T)
  parental <- foreach(x=fns)%do%loadDat(x)
  parental <- rbindlist(parental)
  setkey(parental, id)
  parental[,allele:=gsub("A_Maternal", "A_m", allele)]
  parental[,allele:=gsub("A_Paternal", "A_p", allele)]
  parental[,allele:=gsub("C_Maternal", "C_m", allele)]
  parental[,allele:=gsub("C_Paternal", "C_p", allele)]

  setkey(parental, id)


### iterate through windows to make SNP call per site

  dat <- foreach(i=1:dim(ppl)[1])%dopar%{
    #i<-1
    message(paste(i, dim(ppl)[1], sep=" / "))

    off.tmp <- ppl[i]
    par.temp <- parental[J(off.tmp$V2:off.tmp$V3), nomatch=0]

    A.seg <- par.temp[allele==tstrsplit(off.tmp$V4, "\\|")[[1]]]
    C.seg <- par.temp[allele==tstrsplit(off.tmp$V4, "\\|")[[2]]]

    out.tmp <- data.table(id=A.seg$id, chr=A.seg$chr, pos=A.seg$pos, diplo=off.tmp$V4, clone=off.tmp$clone,
                          geno=paste(A.seg$value, C.seg$value, sep=""))
    return(out.tmp)
  }
  dat <- rbindlist(dat)


### get parental genotypes into it
  parental.ag <- parental[,list(geno=c(paste(value[allele=="A_m"], value[allele=="A_p"], sep=""),
                                       paste(value[allele=="C_m"], value[allele=="C_p"], sep="")),
                                obs.geno=c(A[allele=="A_m"],
                                           C[allele=="C_m"]),
                 clone=c("A", "C"),
                  diplo=c("A_m|A_p", "C_m|C_p")),
            list(chr, pos, id)]

### combine offspring + parental
  #dat <- datp[,-"diplo",with=F]
  datp <- rbind(dat, parental.ag, fill=T, all=T)
  save(datp, file="~/datp.Rdata")
  write.csv(datp, file= "/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all.vitPath.csv")

### load
  datp <- fread(file= "/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all.vitPath.csv")

  setnames(datp, "geno", "phase.geno")



### add in imputed genotypes from Rabbit + original genotype call
  setwd("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm")

  loadDat.io <- function(fn) {
    print(fn)
    #fn <- fns[1]

    ### imputed genotype
      fn <- gsub(".post.csv", "_ImputedGenotype.csv", fn)
      imputed <- fread(fn, skip=1, header=T)
      imputed.l <- melt(imputed, id.vars="marker")

      setnames(imputed.l, c("variable", "marker", "value"), c("id", "clone", "imputedGeno"))
      imputed.l[,id:=as.numeric(as.character(id))]

      setkey(imputed.l, id, clone)
      setkey(datp, id, clone)

      phase.impute <- merge(datp, imputed.l)

    ### get original genotypes
      seqSetFilter(genofile,
                  sample.id=unique(phase.impute$clone),
                  variant.id=unique(phase.impute$id))

      obs.geno <- seqGetData(genofile, "$dosage") ### dosage of the ref allele

      obs.geno.l <- data.table(obs.geno=expand.grid(obs.geno)$Var1,
                              clone=rep(seqGetData(genofile, "sample.id"), dim(obs.geno)[2]),
                              id=rep(seqGetData(genofile, "variant.id"), each=dim(obs.geno)[1]))

      phase.impute.obs <- merge(phase.impute, obs.geno.l, all=T)
      setnames(phase.impute.obs, c("obs.geno.x", "obs.geno.y"), c("consensusParent", c("obs.dosage")))
      return(phase.impute.obs)
    }

    fns <- system("ls /scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/*/*out.post.csv", intern=T)
    pio <- foreach(x=fns)%dopar%loadDat.io(x)
    pio <- rbindlist(pio)

    setnames(pio, "geno", "phase.geno")


    table(pio$phase.geno, pio$imputedGeno)

    table(pio$phase.geno, pio$consensusParent)

    table(pio$phase.geno, pio$obs.dosage)

    table(pio$imputedGeno, pio$consensusParent)

    save(pio, file="~/pio.Rdata")
    write.csv(pio, file= "/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all.vitPath.pio.csv")















### some stuff with teh vit segments
  ### basic summary stats
    ppl.ag <- ppl[,list(nRecomb=length(V4)-1), list(chr, clone)]
    ppl.ag.ag <- ppl.ag[,list(mu=mean(nRecomb)), list(chr)]
    mean(ppl.ag.ag$mu)
    ppl[,A:=tstrsplit(V4, "\\|")[[1]]]
    ppl[,C:=tstrsplit(V4, "\\|")[[2]]]
    table(ppl$A, ppl$chr)
    table(ppl$C)

    table(ppl$A, ppl$C)

  ### distances of paths
    setnames(ppl, "V2", "id")
    setkey(ppl, id)
    setkey(snp.dt, id)
    ppl.start <- merge(ppl, snp.dt)
    setnames(ppl.start, "pos", "start.pos")

    setnames(ppl.start, "id", "V2")
    setnames(ppl.start, "V3", "id")
    setkey(ppl.start, id)
    setkey(snp.dt, id)
    ppl.start.stop <- merge(ppl.start, snp.dt)
    setnames(ppl.start.stop, "pos", "stop.pos")

    pplss <- ppl.start.stop
    pplss[,dist:=stop.pos - start.pos]


  setkey(pplss, clone)
  setkey(sc, clone)
  pplss <- merge(pplss, sc)

  save(pplss, file="~/pplss.Rdata")


  table(pio[chr=="Scaffold_9200_HRSCAF_10757"][pos>8e6]$diplo)
  table(pio[chr=="Scaffold_9200_HRSCAF_10757"][pos>8e6][clone%in%c("A", "C")]$phase.geno,
        pio[chr=="Scaffold_9200_HRSCAF_10757"][pos>8e6][clone%in%c("A", "C")]$clone)

  pio[,list()]
# scp aob2x@rivanna.hpc.virginia.edu:~/pplss.Rdata ~/pplss.Rdata

library(data.table)
library(ggplot2)
load("~/pplss.Rdata")

chromPaint <- ggplot(data=pplss[abs(dist)>100000]) +
geom_segment(aes(y=clone, yend=clone, x=start.pos, xend=stop.pos, color=V4), size=4) +
facet_grid(SC~chr.x, scales="free", space="free_y", shrink=T)

ggsave(chromPaint, file="~/chromPaint.pdf")



  ppl[,method:="consensus.dosage"]
  save(ppl, file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/ppl.consensus.dosage.Rdata")




###
  rm(list=ls())

  load(file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/ppl.consensus.dosage.Rdata")
  ppl.con.d <- ppl

  load(file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/ppl.maxRD.readCounts.Rdata")
  ppl.rc <- ppl

  load(file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/ppl.maxRD.dosage.Rdata")
  ppl[, method:="maxRD.dosage"]
  ppl.d <- ppl

  ppl <- rbindlist(list(ppl.con.d, ppl.rc, ppl.d))

  save(ppl, file="~/3method.ppl.Rdata")

  scp aob2x@rivanna.hpc.virginia.edu:~/3method.ppl.Rdata ~/.

  library(data.table)
  library(ggplot2)
  install.pack
  load(file="~/3method.ppl.Rdata")
  ppl.ag <- ppl[,list(nRecomb=length(V4)), list(chr, clone, method)]

  ggplot(ppl.ag, aes(x=method, y=nRecomb, color=grepl("AxB", clone))) + geom_jitter() + facet_wrap(~chr)

  ppl[,A:=tstrsplit(V4, "\\|")[[1]]]
  ppl[,C:=tstrsplit(V4, "\\|")[[2]]]

  summary(aov(nRecomb~method+chr+clone, ppl.ag))









  #save(ppl, file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/combined_rabbitOut_all_AxC.Rdata")



  setkey(ppl, id)
  setkey(snp.dt, id)

  m <- merge(ppl, snp.dt)

  table(m$chr.x==m$chr.y)

  write.csv(m, file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all_AxC.csv", quote=F, row.names=F)

  write.csv(m, file="/project/berglandlab/alan/all_AxC.csv", quote=F, row.names=F)

















        setnames(pp, names(pp), as.character(pp[1,]))
        pp <- pp[grepl("genotype", marker)]

        ppl <- melt(pp[-(1:3)], id.vars="marker")
        ppl[value!="1",value2:=as.numeric(paste(value, "0", sep=""))]
        ppl[value=="1", value2:=1]

        ppl[,genotype:=last(tstrsplit(marker, "_"))]
        ppl[,clone:=gsub("_genotype[0-9]{1,}", "", marker)]


        ### get only the possible states
          ppl.ag <- ppl[,list(n=sum(value2>.95)), genotype]

          setkey(ppl, genotype)
          setkey(ppl.ag, genotype)

          ppl <- ppl[J(ppl.ag[n>10]$genotype)]

        ### gneotype code
        gc <- fread(fn, skip="Genotype,Code,founder", nrows=10)
        setnames(gc, "Genotype", "genotype")

        ppl <- merge(ppl, gc)
        codes <- data.table(geno=LETTERS[1:4], founder=names(table(ppl$founder))[c(1,2,3,4)])

        ppl <- merge(ppl, codes, by="founder")
        ppl[,chr:=last(tstrsplit(gsub(".all.out.post.csv", "", fn), "/"))]

        setnames(ppl, "variable", "id")
        #setkey(ppl, id)
        #m <- merge(ppl, snp.dt)

        ppl.ag <- ppl[,list(founder=geno[which.max(value2)]), list(id, clone, chr=chr)]
        ppl.ag[,id:=as.numeric(as.character(id))]

      ### check
        chr.i <- "Scaffold_9199_HRSCAF_10755"
        pos.i <- 6229430
        id.i <- snp.dt[chr==chr.i][pos==pos.i]$id
        #ppl.ag[id==id.i]


      ### imputed genotype
        fn <- gsub(".post.csv", "_ImputedGenotype.csv", fn)
        imputed <- fread(fn, skip=1, header=T)
        imputed.l <- melt(imputed, id.vars="marker")

        setnames(imputed.l, c("variable", "marker", "value"), c("id", "clone", "imputedGeno"))
        imputed.l[,id:=as.numeric(as.character(id))]

        setkey(imputed.l, id, clone)
        setkey(ppl.ag, id, clone)

        phase.impute <- merge(ppl.ag, imputed.l)

      ### get original genotypes
        seqSetFilter(genofile,
                    sample.id=unique(phase.impute$clone),
                    variant.id=unique(phase.impute$id))

        obs.geno <- seqGetData(genofile, "$dosage") ### dosage of the ref allele

        obs.geno.l <- data.table(obs.geno=expand.grid(obs.geno)$Var1,
                                clone=rep(seqGetData(genofile, "sample.id"), dim(obs.geno)[2]),
                                id=rep(seqGetData(genofile, "variant.id"), each=dim(obs.geno)[1]))

        phase.impute.obs <- merge(phase.impute, obs.geno.l)
        table(phase.impute.obs[id==id.i]$obs.geno,
              phase.impute.obs[id==id.i]$value)


        phase.impute.obs[id==id.i]























### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### opebn genofiles
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds")

### load hap data
  haps.files <- list.files(path="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/", pattern="*haps")

  haps <- foreach(fi=haps.files, .combine="rbind")%do%{
    #fi <- haps.files[1]
    message(fi)
    tmp <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/", fi, sep="/"))
    tmp[,cloneid:=tstrsplit(V1, " ")[[1]]]
    tmp[,chr:=tstrsplit(V1, " ")[[2]]]
    tmp[,allele1:=tstrsplit(V4, "\\|")[[1]]]
    tmp[,allele2:=tstrsplit(V4, "\\|")[[2]]]

    tmp
  }
  haps[V4=="A_m|C_m", encode:=1]
  haps[V4=="A_p|C_m", encode:=2]
  haps[V4=="A_m|C_p", encode:=3]
  haps[V4=="A_p|C_p", encode:=4]









  seqSetFilter(genofile, unique(haps$chr))

  dt <- data.table(V2=seqGetData(genofile, "variant.id"), startPos=seqGetData(genofile, "position"))
  haps <- merge(haps, dt, by="V2")

  haps[,allele1:=V4]
  haps[,allele2:=V5]

  haps[V4=="A1", allele1:="A"]
  haps[V4=="A2", allele1:="B"]
  haps[V4=="C2", allele1:="D"]
  haps[V4=="C1", allele1:="C"]

  haps[V5=="A1", allele2:="A"]
  haps[V5=="A2", allele2:="B"]
  haps[V5=="C2", allele2:="D"]
  haps[V5=="C1", allele2:="C"]

  haps[,diplo:=apply(cbind(haps$allele1, haps$allele2), 1, function(x) paste(sort(x), collapse=""))]

  prop.table(table(haps$V1, haps$diplo))



  sum(haps[is.na(encode)]$stopPos - haps[is.na(encode)]$startPos)
  sum(haps[!is.na(encode)]$stopPos - haps[!is.na(encode)]$startPos)
