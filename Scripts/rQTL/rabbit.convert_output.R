#module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(SeqArray)
  library(doMC)
  registerDoMC(10)


### make SNP table
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile))
  setkey(snp.dt, id)

### load posterior probabilities
  setwd("/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm")

  loadDat <- function(fn) {
    print(fn)
    #fn <- fns[10]

    ### most likely genotype
      pp <- fread(fn, skip="magicReconstruct-Summary,genoprob,Conditonal genotype probability", header=T, fill=T)

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


      #phase.impute.obs[id==id.i]


  }

  fns <- system("ls /scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/*/*out.post.csv", intern=T)
  ppl <- foreach(x=fns)%dopar%loadDat(x)
  ppl <- rbindlist(ppl)

  #save(ppl, file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/combined_rabbitOut_all_AxC.Rdata")



  setkey(ppl, id)
  setkey(snp.dt, id)

  m <- merge(ppl, snp.dt)

  table(m$chr.x==m$chr.y)

  write.csv(m, file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all_AxC.csv", quote=F, row.names=F)

  write.csv(m, file="/project/berglandlab/alan/all_AxC.csv", quote=F, row.names=F)



























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
