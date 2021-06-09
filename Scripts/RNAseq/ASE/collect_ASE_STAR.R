#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### library
  library(data.table)
  library(foreach)
  library(SeqArray)
  library(tidyverse)

### load GDS file
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

### ase sample table
  sampleTable <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable")
  setnames(sampleTable, "SampleName", "samp")
  sampleTable[,clone:=toupper(clone)]

### superclone file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  sc <- sc[Nonindependent==0 ]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]
  sc[,pond:=toupper(population)]

### replace clone name
  for(i in 1:8) {
    sampleTable[i,clone:=sc[grepl(sampleTable[i]$clone, clone)]$clone]
  }
  sampleTable

### load ASE
  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/ase/", "rnaseq_asereadcounter.star.delim", full.names=T)

  ase <- foreach(i=fn)%do%{
    message(i)
    #i<-fn[1]
    ase.tmp <- fread(i)
    message(dim(ase.tmp))
    ase.tmp[,samp:=tstrsplit(i, "/")%>% last %>% gsub("_rnaseq_asereadcounter.star.delim", "", .)]
    setnames(ase.tmp, c("contig", "position"), c("chr", "pos"))

    ase.tmp[,c("chr", "pos", "refCount", "altCount", "totalCount", "samp"), with=F]
  }
  ase <- rbindlist(ase)

  ase <- merge(ase, sampleTable[,-"id",with=F], by="samp")

### get genotypes
  ase.uniq <- ase[,.N, list(chr, pos)]
  seqSetFilter(genofile, sample.id=unique(sampleTable$clone))
  seqSetFilterPos(genofile, chr=ase.uniq$chr, pos=ase.uniq$pos, intersect=T)

  dosage <- seqGetData(genofile, "$dosage")
  dosage <- as.data.table(t(dosage))
  setnames(dosage, seqGetData(genofile, "sample.id"))
  dosage[,chr:=seqGetData(genofile, "chromosome")]
  dosage[,pos:=seqGetData(genofile, "position")]
  dosage[,id:=seqGetData(genofile, "variant.id")]

  ### check
    prop.table(table(dosage$April_2017_D8_179==dosage$April_2017_D8_349))
    prop.table(table(dosage$April_2017_D8_179==dosage$May_2017_D8_515))

  dosage.long <- melt(dosage, id.vars=c("chr", "pos", "id"),
                      variable.name="clone",
                      value.name="ref_dosage")

### combine ASE and dosage
  setkey(dosage.long, chr, pos, clone)
  setkey(ase, chr, pos, clone)

  ase.geno <- merge(ase, dosage.long, all=T)

  table(ase.geno$altCount>0 & ase.geno$refCount>0, ase.geno$ref_dosage==1)
  table(ase.geno$altCount>0 & ase.geno$refCount==0, ase.geno$ref_dosage==0)
  table(ase.geno$altCount==0 & ase.geno$refCount>0, ase.geno$ref_dosage==2)

### add in annotations
  ase.geno.uniq <- ase.geno[,.N,list(chr, pos, id)]
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=ase.geno.uniq$id)

  tmp <- seqGetData(genofile, "annotation/info/ANN")
  len1 <- tmp$length
  len2 <- tmp$data

  snp.dt1 <- data.table(len=rep(len1, times=len1),
                        ann=len2,
                        id=rep(ase.geno.uniq$id, times=len1))

  snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
  snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]
  snp.dt1.an <- snp.dt1[,list(n=length(class),
                              col= paste(class, collapse=","),
                              class=class[1],
                              genes=gene[1]),
                        list(id=id)]

  snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]

### merge
  ase.geno <- merge(ase.geno, snp.dt1.an, by="id", all=T)
  ase.geno[genes=="Daphnia00796"]
  ase.geno[genes=="Daphnia00787"][!is.na(ref_dosage)][,list(ref=mean(refCount/totalCount, na.rm=T),
                                                            tc=mean(totalCount, na.rm=T), .N), list(clone, ref_dosage)]


### add in phased data

  seqClose(genofile)
  genofile <- seqOpen("/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.gds")

  phased.samps <- seqGetData(genofile, "sample.id")
  setkey(sc, clone)
  phased.samps <- data.table(rClone=c(sc[J(phased.samps)][SC=="A"]$clone, sc[J(phased.samps)][SC=="C"]$clone), SC=c("A", "C"))

  sites.uniq <- ase.geno[,list(.N), list(chr, pos)]

  seqResetFilter(genofile)
  seqSetFilterPos(genofile, chr=sites.uniq$chr, pos=sites.uniq$pos, intersect=T)
  snps.phase.dt <- data.table(chr=seqGetData(genofile, "chromosome"), pos=seqGetData(genofile, "position"), id=seqGetData(genofile, "variant.id"))
  snps.phase.dt <- merge(snps.phase.dt, sites.uniq)
  seqResetFilter(genofile)




  seqSetFilter(genofile,  sample.id=phased.samps$rClone, variant.id=snps.phase.dt$id)
  tmp <- seqGetData(genofile, "genotype")
  o <- foreach(i=c(1,2))%do%{
    tmp.tmp <- as.data.table(t(tmp[i,,]))
    setnames(tmp.tmp, names(tmp.tmp), seqGetData(genofile, "sample.id"))
    tmp.tmp <- cbind(tmp.tmp, snps.phase.dt)
    o <- melt(tmp.tmp, id.vars=c("chr", "pos", "id"), variable.name="clone", value.name="allele")
    setkey(o, chr, pos, id, clone)
    o
  }
  om <- do.call("merge", o)
  om <- merge(om, phased.samps, by.x="clone", by.y="rClone")
  setnames(om, "clone", "rClone")
  setnames(om, "SC", "superclone")

  setkey(ase.geno, chr, pos, superclone)
  setkey(om, chr, pos, superclone)

  ase.geno.phase <- merge(ase.geno, om)
  ase.geno.phase[pos==5204042]

  ase.geno.phase[,xCount:=NA]
  ase.geno.phase[allele.x==1, xCount:=altCount]
  ase.geno.phase[allele.x==0, xCount:=refCount]
  ase.geno.phase[pos==5204042]
  ase.geno.phase[pos==40798]

  ase.geno.phase[genes=="Daphnia00787"][,list(nHet=sum(allele.x!=allele.y), mux=mean(allele.x), muy=mean(allele.y)), list(clone, superclone)]

  table(ase.geno.phase$xCount==ase.geno.phase$altCount, ase.geno.phase$allele.x)


### save
  save(ase.geno.phase, file="~/ase_geno_phase.star.Rdata")

  ### scp aob2x@rivanna.hpc.virginia.edu:~/ase_geno.star.Rdata ~/.

    library(data.table)
    library(ggplot2)

    load("~/ase_geno.star.Rdata")

    ggplot(data=ase.gene[ref_dosage==1][totalCount>500][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")],
      aes(y=refCount/totalCount, x=as.factor(ref_dosage))) + geom_violin() + facet_grid(SC~clone)

    ggplot(data=ase.geno[totalCount>500][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")],
          aes(refCount/totalCount, group=samp, color=samp, fill=samp)) +
    geom_histogram() + facet_wrap(~ref_dosage, scales="free_y")

    ggplot(data=ase.gene[class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant")], aes(y=log10(totalCount), x=class)) +
    geom_boxplot()


    ggplot(data=ase.gene[genes=="Daphnia00796"][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")],
          aes(y=totalCount, x=SC, group=clone, color=clone, fill=clone)) +
    geom_point()



      ggplot(data=ase.geno[genes=="Daphnia00796"][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")],
            aes(refCount/totalCount, group=samp, color=samp, fill=samp)) +
      geom_histogram()

























  ### open genofile
    snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                         pos=seqGetData(genofile, "position"),
                         id=seqGetData(genofile, "variant.id"),
                         numAlleles=seqNumAllele(genofile),
                         ref=seqGetData(genofile, "$ref"),
                         alt=seqGetData(genofile, "$alt"))
    setkey(snp.dt, chr, pos)

  ### clone metadata file
    sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

    sc <- sc[Nonindependent==0 ]

    sc[,SC.uniq:=SC]
    sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]
    sc[,pond:=toupper(population)]

    ### hard filtering of SC
      sc.ag <- sc[LabGenerated==F & Species=="pulex", list(clone=clone[which.max(medrd)][1]), list(SC.uniq, Species)]

  ## iterate through A & C
    geno <- foreach(SC.i=c("A", "C"))%do%{
      seqSetFilter(genofile, sample.id=sc.ag[SC.uniq==SC]$clone)

      tmp <- seqGetData(genofile, "$dosage")

      tmp <- cbind(snp.dt, data.table(ref_dosage=t(tmp)[,1]))
      tmp[,SC:=SC.i]
      tmp
    }
    geno <- rbindlist(geno)
    seqClose(genofile)

### load ASE
  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/ase/", "rnaseq_asereadcounter.delim", full.names=T)

  ase <- foreach(i=fn)%do%{
    #i<-fn[1]
    ase.tmp <- fread(i)
    ase.tmp[,samp:=tstrsplit(i, "/")%>% last %>% gsub("_rnaseq_asereadcounter.delim", "", .)]
    setnames(ase.tmp, c("contig", "position"), c("chr", "pos"))

    ase.tmp[,c("chr", "pos", "refCount", "altCount", "totalCount", "samp"), with=F]
  }
  ase <- rbindlist(ase)

  sampleTable <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable")
  setnames(sampleTable, "SampleName", "samp")
  ase <- merge(ase, sampleTable[,-"id",with=F], by="samp")
  setnames(ase, "superclone", "SC")
  setkey(ase, chr, pos, SC)
  setkey(geno, chr, pos, SC)

  ase <- merge(ase, geno)

### get annotations
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

  ### load meta-data file
    samps <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  ### load GDS file
    genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

  ### filter
    ase.uniq <- ase[,list(.N), list(chr, pos, id)]
    seqSetFilterPos(genofile, chr=ase.uniq$chr, pos=ase.uniq$pos)

    tmp <- seqGetData(genofile, "annotation/info/ANN")
    len1 <- tmp$length
    len2 <- tmp$data

    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          id=rep(ase.uniq$id, times=len1))

    snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
    snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]
    snp.dt1.an <- snp.dt1[,list(n=length(class),
                                col= paste(class, collapse=","),
                                class=class[1],
                                genes=gene[1]),
                          list(id=id)]

    snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]

  ### tack in
    ase.gene <- merge(ase, snp.dt1.an, by="id")

### save
  save(ase.gene, file="~/ase_gene.Rdata")

### scp aob2x@rivanna.hpc.virginia.edu:~/ase_gene.Rdata ~/.

  library(data.table)
  library(ggplot2)

  load("~/ase_gene.Rdata")

  ggplot(data=ase.gene[ref_dosage==1][totalCount>500][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")],
    aes(y=refCount/totalCount, x=as.factor(ref_dosage))) + geom_violin() + facet_grid(SC~clone)

  ggplot(data=ase.gene[totalCount>500][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")],
        aes(refCount/totalCount, group=samp, color=samp, fill=samp)) +
  geom_histogram() + facet_wrap(~ref_dosage, scales="free_y")

  ggplot(data=ase.gene[class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant")], aes(y=log10(totalCount), x=class)) +
  geom_boxplot()


  ggplot(data=ase.gene[genes=="Daphnia00796"][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")],
        aes(y=totalCount, x=SC, group=clone, color=clone, fill=clone)) +
  geom_point()



    ggplot(data=ase.gene[genes=="Daphnia00796"][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")],
          aes(refCount/totalCount, group=samp, color=samp, fill=samp)) +
    geom_histogram()



    snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                         pos=seqGetData(genofile, "position"),
                         id=seqGetData(genofile, "variant.id"),
                         numAlleles=seqNumAllele(genofile),
                         key="chr")
    setkey(snpFilter, chr, pos)
    setkey(snp.dt, chr, pos)

    snp.dt <- merge(snpFilter, snp.dt)
