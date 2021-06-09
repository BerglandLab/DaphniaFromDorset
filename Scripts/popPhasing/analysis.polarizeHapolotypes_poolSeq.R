bcftools view \
-O v \
-o /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf \
/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.bcf

#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(foreach)
  library(SeqArray)
  library(SeqVarTools)

### convert to GDS
  vcf.fn <- "/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf"
  gds.fn <- "/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.gds"

  #seqVCF2GDS(vcf.fn, gds.fn)

### open genofile
  genofile <- seqOpen(gds.fn)
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

  ### sc per year
    sc.peryear <- sc[,list(.N), list(year, SC.uniq, pond)]

  ### hard filtering of SC
    sc.ag <- sc[LabGenerated==F & Species=="pulex", list(clone=clone[which.max(medrd)][1]), list(SC.uniq, Species)]
    #sc.ag[,pond:=toupper(pond)]

### load poolseq data
  load("/project/berglandlab/alan/gprime_peaks.replicates.250K.05.Rdata")
  setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))
  setkey(gprime, chr, pos,rep)

### function
  polarize <- function(i, window) {
    #window<- 10000; i<-10
    chr.i<-peaks[i]$CHROM; start<-peaks[i]$posMaxGprime-window; stop<-peaks[i]$posMaxGprime+window; rep.i=peaks[i]$rep

    pool.tmp <- gprime[J(data.table(chr=chr.i, pos=start:stop, rep=rep.i, key="chr,pos,rep")), nomatch=0]
    setkey(pool.tmp, chr, pos)
    setkey(snp.dt, chr, pos)
    unique(snp.dt[pool.tmp]$id)
    seqSetFilter(genofile, variant.id=unique(snp.dt[pool.tmp]$id))

    tmp <- as.data.table(getGenotype(genofile))
    tmp[,sample.id:=seqGetData(genofile, "sample.id")]

    dat.phase <- melt(tmp, id.vars="sample.id", variable.name="variant.id", value.name="geno")
    dat.phase[,allele1:=tstrsplit(geno, "\\|")[[1]]]
    dat.phase[,allele2:=tstrsplit(geno, "\\|")[[2]]]

    # dat.phase[sample.id=="April_2017_D8_213"] ## A
    # dat.phase[sample.id=="April_2017_D8_151"] ## C

    dat.phase <- melt(dat.phase[,c("sample.id", "variant.id", "allele1", "allele2")],
                    id.vars=c("sample.id", "variant.id"),
                    value.vars=c("allele1", "allele2"),
                    value.name="allele")

    dat.phase[,variable:=gsub("allele", "", variable)]
    dat.phase[,haplo:=paste(sample.id, variable, sep=".")]
    setnames(dat.phase, "variant.id", "id")
    dat.phase[,id:=as.numeric(as.character(id))]

    setkey(dat.phase, id)
    setkey(snp.dt, id)


    dat.phase <- merge(dat.phase, snp.dt)
    dat.phase[,rep:=rep.i]

    setkey(dat.phase, chr, pos, rep)
    setkey(gprime, chr, pos, rep)

    m <- merge(dat.phase, gprime)

    m[sign(deltaSNP)== -1 &  allele==0, concord:="male"]
    m[sign(deltaSNP)== -1 &  allele==1, concord:="pe"]
    m[sign(deltaSNP)==  1 &  allele==0, concord:="pe"]
    m[sign(deltaSNP)==  1 &  allele==1, concord:="male"]
    m <- merge(m, sc.ag, by.x="sample.id", by.y="clone")

    # dcast(m[SC.uniq=="C"], chr+pos+SC.uniq~variable, value.var="concord")
    # dcast(m[SC.uniq=="A"], chr+pos~variable, value.var="concord")

    m.ag <- dcast(m, chr+pos+SC.uniq~variable, value.var="concord")
    setnames(m.ag, c("1","2"), c("allele1", "allele2"))
    m.ag.ag <- m.ag[,list(n.male_male=sum(allele1=="male" & allele2=="male")/length(allele1),
                  n.male_pe=sum((allele1=="male" & allele2=="pe") | (allele2=="male" & allele1=="pe"))/length(allele1),
                  n.pe_pe=sum(allele1=="pe" & allele2=="pe")/length(allele1)),
              list(SC.uniq)]
    setkey(m.ag.ag, SC.uniq)
    setkey(sc.ag, SC.uniq)

    m.ag.ag <- merge(m.ag.ag, sc.ag)
    m.ag.ag[,pond:=toupper(tstrsplit(clone, "_")[[3]])]
    m.ag.ag[,qtl:=i]
    m.ag.ag

  }


  qtl.polar <- foreach(i=1:14, .combine="rbind")%do%polarize(i, window=20000)

  setkey(qtl.polar, SC.uniq)
  setkey(sc.peryear, SC.uniq)
  qtl.polar <- merge(qtl.polar, sc.peryear)

  save(qtl.polar, file="~/qtl_polar.Rdata")

### download and plot
  scp aob2x@rivanna.hpc.virginia.edu:~/qtl_polar.Rdata ~/.

library(data.table)
library(ggplot2)
library(patchwork)
  load("~/qtl_polar.Rdata")


#### abundance genotype
  ab <- qtl.polar[,list(geno=c("male_male", "male_pe", "pe_pe")[which.max(c(n.male_male, n.male_pe, n.pe_pe))],
                        size=N[1]), list(clone, pond=pond.y, qtl, sc=SC.uniq, year=year)]

  abr <- ab[,list(geno=rep(geno, size)), list(clone, pond, sc, qtl, year)]

  abrf <- abr[,list(male_freq=(2*sum(geno=="male_male") + sum(geno=="male_pe"))/(2*length(geno)), n=2*length(geno)), list(pond, year, qtl)]
  abrf[,se:=male_freq*(1-male_freq)/sqrt(n)]
  #abrf[,pond:=factor(pond, levels=c("DBUNK", "D8", "DCAT"))]

  good <- ggplot(data=abrf[pond%in%c("D8", "DBUNK", "DCAT")][year>2016], aes(x=as.factor(year), y=male_freq, group=pond, color=pond)) +
  geom_errorbar(aes(ymin=male_freq-2*se, ymax=male_freq+2*se)) +
  geom_point() + geom_line() + facet_grid(~qtl) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  good / bad
  table(abr[qtl==2]$geno, abr[qtl==2]$pond)
  table(big=ab[pond=="D8"][qtl==8]$size>=2, ab[pond=="D8"][qtl==8]$geno)


#### some pond level stuff
  qtl.polar[,list()]
  qtl.polar[qtl==10][pond=="D8"]
  qtl.polar[qtl==8][pond=="DBUNK"]

  qtl.polar[qtl==8][]
  ### ponds
  qtl.polar.ag <- qtl.polar[,list(male_male=sum(n.male_male>.75), male_pe=sum(n.male_pe>.75), pe_pe=sum(n.pe_pe>.75)), list(pond, qtl)]

  qtl.polar.ag[,male_male_frac:=male_male/(male_male + male_pe + pe_pe)]
  qtl.polar.ag[,male_pe_frac:=  male_pe  /(male_male + male_pe + pe_pe)]
  qtl.polar.ag[,pe_pe_frac:=    pe_pe    /(male_male + male_pe + pe_pe)]

  qtl.polar.ag.l <- melt(qtl.polar, id.vars=c("SC.uniq", "clone", "pond", "qtl", "Species"), variable.name="geno", value.name="pr")
  qtl.polar.ag.l.ag <- qtl.polar.ag.l[,list(geno=geno[which.max(pr)]), list(SC.uniq, clone, pond, qtl)]

  qtl.polar.ag.l.ag.ag <- qtl.polar.ag.l.ag[,list(n=(2*sum(geno=="n.male_male") + 1*sum(geno=="n.male_pe"))/(2*length(geno))),
                                             list(pond, qtl)]

   sum <- qtl.polar.ag.l.ag[,list(prop=c(sum(geno=="n.male_male") / length(geno),
                                         sum(geno=="n.male_pe") / length(geno),
                                         sum(geno=="n.pe_pe") / length(geno)),
                                  n=c(sum(geno=="n.male_male"),
                                      sum(geno=="n.male_pe"),
                                      sum(geno=="n.pe_pe")),
                                  N=length(geno),
                                  geno=c("male_male", "male_pe", "pe_pe")),
                      list(pond, qtl)]

    sum[qtl==10][pond=="DBUNK"]



    sum.ag <- sum[,list(diff=prop[pond=="DBUNK"]-prop[pond=="D8"],
              mean=(prop[pond=="DBUNK"]+prop[pond=="D8"])/2), list(qtl, geno)]
    sum.ag[,cv:=diff/mean]





    ggplot(data=sum.ag, aes(x=geno, y=cv)) + geom_point() + facet_grid(~qtl)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





ggplot(data=sum[pond%in%c("DBUNK")], aes(x=geno, y=prop, group=pond, color=pond)) +
geom_line() + facet_grid(~qtl) +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


 male <- ggplot(qtl.polar.ag.l.ag.ag[SC.uniq%in%c("A", "C")], aes(x=pond, y=propMaleAlleles, group=qtl)) +
         #geom_point(aes(color=pond)) +
         geom_line() + facet_grid(~qtl, scales="free")



  male <- ggplot(qtl.polar.ag.l.ag.ag[pond%in%c("D8", "DBUNK")], aes(x=pond, y=propMaleAlleles, group=qtl)) +
          #geom_point(aes(color=pond)) +
          geom_line() + facet_grid(~qtl)











  sc <- ggplot(qtl.polar.ag.l.ag[grepl("18004", clone) | SC.uniq%in%c("A", "B", "C") | grepl("po", SC.uniq)], aes(x=qtl, y=SC.uniq, fill=geno)) + facet_grid(~qtl, scales="free") + geom_tile()
  pond <- ggplot(qtl.polar.ag.l.ag[pond%in%c("D8", "DBUNK")], aes(x=qtl, y=SC.uniq, fill=geno)) + facet_grid(pond~qtl, scales="free") + geom_tile()

  ggsave(sc / pond , file="~/polarized_haplotypes.pdf", h=10, w=15)


  qtl.polar.ag[qtl==10][pond%in%c("D8", "DBUNK", "DCAT")]
  qtl.polar.ag[qtl==8][pond%in%c("D8", "DBUNK", "DCAT")]

  qtl.polar[grepl("18004", clone) | SC.uniq%in%c("A", "B", "C")][qtl==6]

### sc
  qtl.polar.ag <- qtl.polar[,list(male_male=sum(n.male_male>.75), male_pe=sum(n.male_pe>.75), pe_pe=sum(n.pe_pe>.75)), list(pond, SC.uniq)]

  qtl.polar.ag[,male_male_frac:=male_male/(male_male + male_pe + pe_pe)]
  qtl.polar.ag[,male_pe_frac:=  male_pe  /(male_male + male_pe + pe_pe)]
  qtl.polar.ag[,pe_pe_frac:=    pe_pe    /(male_male + male_pe + pe_pe)]

  qtl.polar.ag[SC.uniq%in%c("A", "B", "C")



  qtl.polar[SC.uniq=="A"]
  qtl.polar[SC.uniq=="C"]

  qtl.polar[,n.male_male:=round(n.male_male, 1)]
  qtl.polar[,n.male_pe:=round(n.male_pe, 1)]
  qtl.polar[,n.pe_pe:=round(n.pe_pe, 1)]
  qtl.polar[grepl("18004", clone)]

  qtl.polar.ag <- qtl.polar[,list(.N), list(frac_male, pond, qtl)][order(pond)]
  qtl.polar.ag[pond%in%c("D8", "DBUNK")][order(qtl)]


  qtl.polar.ag <- qtl.polar[,list(which.max(c(n.male_n.male, n.male_pe, n.pe_pe)))]

  ggplot(data=qtl.polar.ag[pond%in%c("D8", "DBUNK")], aes(x=frac_male*2, y=N, group=pond, color=pond)) +
  geom_line() + facet_wrap(~qtl) + geom_point()



  qtl.polar.ag <- qtl.polar[,list(.N), list(frac_male, SC.uniq, qtl)]
  qtl.polar.ag[SC.uniq%in%c("A", "B", "C")][order(qtl)]

  ggplot(data=qtl.polar.ag[SC.uniq%in%c("A", "B", "C")], aes(fill=as.factor(round(frac_male*2)), y=SC.uniq, x=qtl)) +
  geom_tile() + facet_grid(~qtl, scales="free")


  qtl.polar.ag <- qtl.polar[frac_male%in%c(0,.5,1)][,list(.N), list(frac_male, pond, qtl)][order(pond)]
qtl.polar.ag[pond%in%c("D8", "DBUNK")][order(qtl)]

qtl.polar.ag <- qtl.polar[frac_male%in%c(0,.5,1)][,list(.N), list(frac_male, SC.uniq, qtl)]
qtl.polar.ag[SC.uniq%in%c("A", "B", "C")][order(qtl)]
