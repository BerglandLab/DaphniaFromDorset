
### copy data
  scp aob2x@rivanna.hpc.virginia.edu:~/cdlo_250K.Rdata ~/.
  scp aob2x@rivanna.hpc.virginia.edu:~/regions.Rdata ~/.

### libraries
  library(ggplot2)
  library(data.table)
  library(ggtree)
  library(ape)

### load data
  load("~/cdlo_250K.Rdata")
  load("~/gprime_peaks.replicates.Rdata")
  load("~/regions.Rdata")

  setnames(peaks, "CHROM", "chr")

  cdl.o <- merge(cdl.o, regions, by="window")
  setnames(cdl.o, "N", "numMissing")
  cdl.o[,fracMissing:=numMissing/size]

  cdl.qtl <- merge(cdl.qtl, regions, by="window")
  setnames(cdl.qtl, "N", "numMissing")

  cdl.qtl[,fracMissing:=numMissing/size]


  ### cdl.o[,fracMissing:=0]; cdl.qtl[,fracMissing:=0]

### genome-wide distributions
  cdl.o[fracMissing<.25,list(n=sum(n)), list(window, sp.group, pond.group)]

  cdl.o.genomeWideDist <- cdl.o[fracMissing<.25,list(n=sum(n)), list(cd_bin, pond.group, sp.group)]

  ggplot(data=cdl.o.genomeWideDist, aes(x=cd_bin, y=log10(n), group=pond.group, color=pond.group)) +
  geom_line() + facet_grid(~sp.group+pond.group, scales="free_y") +
  geom_vline(xintercept=0.004)


### distance probabilities for each QTL:
  cdl.o.distProb <- cdl.o[fracMissing<.75, list(max=max(cd_bin[n>0])), list(window, pond.group, sp.group)]
  cdl.o.distProb[,region:=tstrsplit(window, ":")[[2]]]
  cdl.o.distProb[,start:=as.numeric(tstrsplit(region, "-")[[1]])]
  cdl.o.distProb[,stop:=as.numeric(tstrsplit(region, "-")[[2]])]

  cdl.qtl.ag <- cdl.qtl[group.x%in%c("A", "C") & group.y%in%c("A", "C")][i1!=i2][,list(max=max(cd)), list(window, pond.group, sp.group)]
  cdl.qtl.ag[,region:=tstrsplit(window, ":")[[2]]]
  cdl.qtl.ag[,start:=as.numeric(tstrsplit(region, "-")[[1]])]
  cdl.qtl.ag[,stop:=as.numeric(tstrsplit(region, "-")[[2]])]

  qtl.age <- foreach(i=1:dim(cdl.qtl.ag)[1])%do%{
                cdl.o.distProb[!(start<=cdl.qtl.ag[i]$start & stop>=cdl.qtl.ag[i]$stop),
                           list(pr=mean(max >= cdl.qtl.ag[i]$max, na.rm=T),window=cdl.qtl.ag[i]$window),
                            list(sp.group, pond.group)]

              }
  qtl.age <- rbindlist(qtl.age)

  qtl.age[pond.group=="DWT-DWT"][order(pr)]
  qtl.age[pond.group=="DWT-D10"][order(pr)]
  qtl.age[pond.group=="DWT-W"][order(pr)]

  ggplot(data=qtl.age, aes(x=pond.group, y=pr, group=window, color=window))+ geom_line()

  ggplot(data=cdl.o.distProb[pond.group!="all"], aes(max)) +
  geom_histogram(bins=200) + facet_grid(pond.group~sp.group) +
  geom_point(data=cdl.qtl.ag, aes(x=max, y=5), color="red")

### manhattan plot
  cdl.o[,chr:=tstrsplit(window, ":")[[1]]]
  cdl.o[,range:=tstrsplit(window, ":")[[2]]]
  cdl.o[,start:=as.numeric(tstrsplit(range, "-")[[1]])]
  cdl.o[,stop:=as.numeric(tstrsplit(range, "-")[[2]])]
  cdl.o[,mid:=start/2 + stop/2]


  cdl.o.manhattanPlot <- cdl.o[,list(mu=sum(cd_bin*n, na.rm=T)/sum(n), min=min(cd_bin), max=max(cd_bin)),
                      list(chr, mid, sp.group, pond.group)]


  mp_small <- ggplot(data=cdl.o.manhattanPlot[!is.na(pond.group)], aes(x=mid, y=max, color=chr)) +
  geom_vline(data=peaks, aes(xintercept=posMaxGprime)) +
  geom_line(size=1) +
  facet_grid(sp.group+pond.group~chr, scales="free")

  ggsave(mp_small, file="~/mp_small.png", height=10, w=20)


### phylogenetic tree
  cdl.o.distProb[,chr:=tstrsplit(window, ":")[[1]]]
  cdl.o.distProb[pond.group=="DWT-DWT"][,list(window=window[which.max(max)], max=max(max)), list(chr)

  ### QTL 8
    localPeak <- "Scaffold_2158_HRSCAF_2565:1750002-2000001"
    QTLPeak <- "Scaffold_2158_HRSCAF_2565:1645971-1895970"

  ### QTL10

    ### local peak
      njo <- cdl.tree$"Scaffold_2158_HRSCAF_2565:1750002-2000001"
      njo <- root(njo, "pulicaria.1")

      d <- data.table(label=njo$tip.label)
      d[,pond:=tstrsplit(label, "_")[[3]]]
      d[,A:=ifelse(grepl("April_2017_D8_213", label), "A", "")]
      d[,C:=ifelse(grepl("April_2017_D8_151", label), "C", "")]
      d <- as.data.frame(d)


      p.local <- ggtree(njo) %<+% d +
      geom_tiplab(aes(label=A)) +
      geom_tiplab(aes(label=C)) +
      geom_tiplab(aes(label=pond, color=pond), offset=0, angle=90, align=T)  +
      coord_flip() +
      theme(legend.position = "none") + ggtitle("local peak")
      #p

    ### qtl peak
      njo <- cdl.tree$"Scaffold_2158_HRSCAF_2565:1645971-1895970"
      njo <- root(njo, "pulicaria.1")

      d <- data.table(label=njo$tip.label)
      d[,label:=gsub("Dcat", "DCat", label)]
      d[,pond:=tstrsplit(label, "_")[[3]]]
      d[,A:=ifelse(grepl("April_2017_D8_213", label), "A", "")]
      d[,C:=ifelse(grepl("April_2017_D8_151", label), "C", "")]

      d[grepl("March20_2018_Dcat_2", label)]

      d[A!="" | C!="" | pond=="DCat" | pond=="Dcat", clean:=label]
      d[grepl("March20_2018_DCat_2", label)]

      d[is.na(clean), clean:=""]
      d <- as.data.frame(d)


      p.qtl <- ggtree(njo) %<+% d +
      geom_tiplab(aes(label=A)) +
      geom_tiplab(aes(label=C)) +
      geom_tiplab(aes(label=pond, color=pond), offset=0, angle=90, align=T)  +
      coord_flip() +
      theme(legend.position = "none")+ ggtitle("qtl-peak")
      p.qtl


      ggtree(njo) %<+% d +
      geom_tiplab(aes(label=A)) +
      geom_tiplab(aes(label=C)) +
      geom_tiplab(aes(label=clean, color=pond), offset=0, angle=0, align=F)  +
      geom_tiplab(aes(label=pond, color=pond), offset=0, angle=0, align=T)  +
      theme(legend.position = "none")+ ggtitle("QTL10")


    ### combined
      p.qtl + p.local


    ### cophyplot
      njo.local <- cdl.tree$"Scaffold_2158_HRSCAF_2565:1750002-2000001"
      njo.local <- root(njo.local, "pulicaria.1")

      njo.qtl <- cdl.tree$"Scaffold_2158_HRSCAF_2565:1645971-1895970"
      njo.qtl <- root(njo.qtl, "pulicaria.1")

      mat

      cophyloplot(njo.local, njo.qtl, assoc=cbind(matrix(njo.local$tip.label, ncol=1), matrix(njo.local$tip.label, ncol=1)),
                  show.tip.label = F, length.line = 4, space = 28, gap = 3, type="cladogram")



  ######
  ## QTL 10 ##
      ### qtl peak
      #  njo <- cdl.tree$"Scaffold_2217_HRSCAF_2652:5073222-5323221"
        njo <- cdl.tree$"Scaffold_9197_HRSCAF_10753:2452444-2702443"
        njo <- root(njo, "pulicaria.1")


        d <- data.table(label=njo$tip.label)
        d[,label:=gsub("Dcat", "DCat", label)]
        d[,pond:=tstrsplit(label, "_")[[3]]]


        d[,A:=ifelse(grepl("April_2017_D8_213", label), "A", "")]
        d[,C:=ifelse(grepl("April_2017_D8_151", label), "C", "")]

        d[grepl("March20_2018_Dcat_2", label)]

        d[A!="" | C!="" | pond=="DCat" | pond=="Dcat", clean:=label]
        d[grepl("March20_2018_Dcat_2", label)]

        d[is.na(clean), clean:=""]
        d <- as.data.frame(d)

        ggtree(njo) %<+% d +
        geom_tiplab(aes(label=A)) +
        geom_tiplab(aes(label=C)) +
        geom_tiplab(aes(label=clean, color=pond), offset=0, angle=0, align=F)  +
        geom_tiplab(aes(label=pond, color=pond), offset=0, angle=0, align=T)  +
        theme(legend.position = "none")+ ggtitle("QTL10")





  #mp <- ggplot(data=cdl.o[!is.na(sp.group)][!is.na(pond.group)][order(n)], aes(x=mid, y=cd_bin, color=log10(n))) + geom_point () +
  #facet_grid(sp.group+pond.group~chr, scales="free_x")
#
  #ggsave(mp, file="~/mp.png", height=10, w=20)
#

  cdl.o.ag <- cdl.o[,list(mu=sum(cd_bin*n, na.rm=T)/sum(n), min=min(cd_bin), max=max(cd_bin)),
                      list(chr, mid, sp.group, pond.group)]


  cdl.o.ag <- cdl.o.ag[,list(mu=mu, min=min, max=max, xi=rank(max, ties="random"), chr=chr, mid=mid), list(sp.group, pond.group)]

  cdl.o.ag[,target:=FALSE]
  cdl.o.ag[round(mid-1.5, -2)!=(mid-1.5), target:=TRUE]

  ggplot(data=cdl.o.ag[!is.na(sp.group)][!is.na(pond.group)][order(target)],
        aes(x=xi, y=max, color=as.factor(target))) +
  geom_point() +
  geom_point(data=cdl.o.ag[!is.na(sp.group)][!is.na(pond.group)][target==T][chr=="Scaffold_2217_HRSCAF_2652"],
        aes(x=xi, y=max), color="black") +
  facet_grid(sp.group+pond.group~.)







  mp_small <- ggplot(data=cdl.o.ag[!is.na(pond.group)], aes(x=mid, y=max, color=chr)) +
  geom_vline(data=peaks, aes(xintercept=posMaxGprime)) +
  geom_line(size=1) +
  facet_grid(sp.group+pond.group~chr, scales="free")

  ggsave(mp_small, file="~/mp_small.png", height=10, w=20)

### all QTL
  o.plot <- ggplot() +
  geom_violin(data=cdl.qtl[!is.na(sp.group)][!is.na(pond.group)][sp.group!="pulicaria-obtusa"],
              aes(y=cd, x=sp.group, color=interaction(sp.group, pond.group), fill=pond.group)) +
  geom_point(data=cdl.qtl[group.x%in%c("A", "C") & group.y%in%c("A", "C")][i1!=i2][sp.group!="pulicaria-obtusa"],
            aes(y=cd, x=sp.group, shape=interaction(group.x, group.y)),
            position = position_nudge(x = -0.5),
            size=2) +
  facet_wrap(~window)

  ggsave(o.plot, file="~/qtl_age.pdf", height=20, width=20)


### single QTL
cdl.genome[, list(TT=sum(n[cd_bin>=0.0045432596]),
                  FF=sum(n[cd_bin<0.0045432596]),
                  pr=sum(n[cd_bin>=0.0045432596]) /sum(n)),
             list(sp.group, pond.group)]

 cdl.genome.big <- cdl.genome[, list(cd=rep(cd_bin, n)),
              list(sp.group, pond.group)]

cdl.genome.big[,list(pr=mean(cd>0.007475457)), list(sp.group, pond.group)]


cdl.qtl[group.x%in%c("A", "C") & group.y%in%c("A", "C")][grepl("Scaffold_2217_HRSCAF_2652", window)]



  ggplot() +
  geom_violin(data=cdl.qtl[!is.na(sp.group)][!is.na(pond.group)][grepl("Scaffold_2217_HRSCAF_2652", window)],
              aes(y=cd, x=sp.group,fill=pond.group)) +
  geom_point(data=cdl.qtl[group.x%in%c("A", "C") & group.y%in%c("A", "C")][i1!=i2][grepl("17459", window)],
            aes(y=cd, x=sp.group, shape=interaction(group.x, group.y)),
            position = position_nudge(x = -0.5),
            size=2) +
  #geom_hline(data=cdl.genome.ag, aes(yintercept=cd, group=interaction(sp.group, pond.group), color=interaction(sp.group, pond.group))) +
  facet_wrap(~window)




  ggplot() +
  geom_violin(data=cdl.qtl[!is.na(sp.group)][!is.na(pond.group)][grepl("17459", window)],
              aes(y=cd, x=sp.group,fill=pond.group)) +
  geom_point(data=cdl.qtl[group.x%in%c("A", "C") & group.y%in%c("A", "C")][i1!=i2][grepl("17459", window)],
            aes(y=cd, x=sp.group, shape=interaction(group.x, group.y)),
            position = position_nudge(x = -0.5),
            size=2) +
  geom_boxplot(data=cdl.genome[grepl("17459", window)],
            aes(y=cd_mean, x=sp.group, group=interaction(sp.group, pond.group), color=interaction(sp.group, pond.group))) +
  facet_wrap(~window)






cdl.qtl[group.x%in%c("A", "C") & group.y%in%c("A", "C")][window=="Scaffold_9199_HRSCAF_10755:4964452-5014451"]

  p <- ggtree(njo) %<+% d +
  geom_tiplab(aes(label=A)) +
  geom_tiplab(aes(label=C)) +
  geom_tiplab(aes(label=pond, color=pond), offset=1, angle=90, align=T)  +
  coord_flip() +
  theme(legend.position = "none")
  p




njo <- cdl.tree$"Scaffold_2158_HRSCAF_2565:1645971-1895970"
njo <- root(njo, "pulicaria.1")








### load reference genome
  read.FASTA(file, type = "DNA")




  genofile <- seqOpen("/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.filter.snpsvarpulex.whatshap.gds")

### snp.dt
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile),
                       key="chr")

### clone metadata file
  sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  sc <- sc[Nonindependent==0]

  sc[,SC.uniq:=SC]
  sc[SC=="OO", SC.uniq:=paste(SC, SCnum, sep="")]

### hard filtering of SC
  sc.ag <- sc[Species=="pulex" & LabGenerated==F, list(clone=clone[which.max(medrd)][1]), list(SC.uniq)]

### ij

ij <- foreach(i=c(1:(dim(mat)[1]-1)), .combine="rbind")%do%{
  foreach(j=c((i+1):dim(mat)[1]), .combine="rbind")%do%{
    data.table(i=i, j=j)
  }
}

### funcitons
  getMatn <- function(samples, chr.x, start, stop) {
    ### samples <- sc.ag[SC.uniq%in%c("A", "C")]$clone; chr.x="Scaffold_2217_HRSCAF_2652"; start=(5198221-10000); stop=(5198221+10000)
    ### samples <- sc.ag$clone; chr.x="Scaffold_2217_HRSCAF_2652"; start=(5198221-10000); stop=(5198221+10000)

      ids <- snp.dt[chr==chr.x][pos>=start & pos<=stop]$id

    ### get genotypes & haplotypes
      seqResetFilter(genofile)
      seqSetFilter(genofile, sample.id=samples, variant.id=ids)
      tmp <- as.data.table(getGenotype(genofile))
      tmp[,sample.id:=seqGetData(genofile, "sample.id")]

      dat.phase <- melt(tmp, id.vars="sample.id", variable.name="variant.id", value.name="geno")

      dat.phase[,allele1:=tstrsplit(geno, "\\|")[[1]]]
      dat.phase[,allele2:=tstrsplit(geno, "\\|")[[2]]]

    ### convert back to wide
      mat <- t(as.matrix(dcast(dat.phase, variant.id~sample.id, value.var=c("allele1", "allele2"))))
      dimnames(mat)[[2]] <- mat[1,]
      mat <- (mat[-1,])
      matn <- apply(mat, 2, as.numeric)

      dimnames(matn)[[1]] <- dimnames(mat)[[1]]
      dimnames(matn)[[2]] <- dimnames(mat)[[2]]
      # save(matn, file="~/phased_matrix.Rdata")

    ### return
      return(matn)
  }


### load peaks file
  load("/scratch/aob2x/daphnia_hwe_sims/gprime_peaks.replicates.250K.05.Rdata")

### get mats
  mats <- foreach(i=1:dim(peaks)[1])%do%{
    print(i)
    matn <- getMatn(samples = sc.ag$clone,
                    chr.x=peaks[i]$CHROM,
                    start=(peaks[i]$posMaxGprime-10000),
                    stop=(peaks[i]$posMaxGprime+10000))
  }

### save
  save(mats, file="~/mats.Rdata")



### scp aob2x@rivanna.hpc.virginia.edu:~/mats.Rdata ~/.

















    ### output


      o <- foreach(x=1:dim(ij)[1], .combine="rbind")%do%{
        i <- ij[x]$i
        j <- ij[x]$j
        data.table(haplo1=dimnames(mat)[[1]][i],
                   haplo2=dimnames(mat)[[1]][j],
                   dxy=mean(matn[i,]!=matn[j,]),
                   n=dim(matn)[2],
                   chr=chr.x, start=start, stop=stop)

      }



  }

### sliding window
  size.bp <- 250000
  step.bp <- 50000

  win.df <- foreach(chr.i=unique(snp.dt$chr), .combine="rbind")%do%{
    data.table(chr=chr.i,
                start=seq(from=min(snp.dt[chr==chr.i]$pos),
                          to=max(snp.dt[chr==chr.i]$pos),
                          by=step.bp),
                stop=seq(from=min(snp.dt[chr==chr.i]$pos),
                          to=max(snp.dt[chr==chr.i]$pos),
                          by=step.bp) + size.bp)
  }

### iterate through
  o <-foreach(i=1:dim(win.df)[1], .errorhandling="remove")%do%{
    message(paste(i, dim(win.df)[1], sep=" / "))
    #i<-100
#    getDists(samples=sc.ag[SC.uniq%in%c("A", "C")]$clone, chr.x=win.df[i]$chr, start=win.df[i]$start, stop=win.df[i]$stop)
    getDists(samples=sc.ag$clone, chr.x=win.df[i]$chr, start=win.df[i]$start, stop=win.df[i]$stop)

  }
  o <- rbindlist(o)

### save
  save(o, file="~/AC_phased_pwd.Rdata")



### scp aob2x@rivanna.hpc.virginia.edu:~/AC_phased_pwd.Rdata ~/.
scp aob2x@rivanna.hpc.virginia.edu:~/phased_matrix.Rdata ~/.


library(data.table)
library(ggplot2)

load("~/AC_phased_pwd.Rdata")

ggplot() +
geom_point(data=o, aes(x=n, y=dxy)) +
geom_point(data=o[chr=="Scaffold_2217_HRSCAF_2652"][start<=5198221 & stop>=5198221], aes(x=n, y=dxy), color="red") +
facet_grid(~haplo1+haplo2)


ggplot() +
geom_point(data=o, aes(x=n, y=dxy)) +
geom_point(data=o[chr=="Scaffold_7757_HRSCAF_8726"][start<=8441175 & stop>=8441175], aes(x=n, y=dxy), color="red") +
facet_grid(~haplo1+haplo2)


library(data.table)
library(ape)
library(ggtree)

load("~/phased_matrix.Rdata")



dm <- dist(matn, diag=T, upper=T)
njo <- nj(dm)
njo <- root(njo, outgroup="allele1_Spring_2016_W6_6.3")
d <- data.table(label=njo$tip.label)
d[,pond:=tstrsplit(label, "_")[[4]]]
d[,A:=ifelse(grepl("April_2017_D8_213", label), "A", "")]
d[,C:=ifelse(grepl("April_2017_D8_151", label), "C", "")]
d <- as.data.frame(d)

ggplot(njo) + geom_tree() +
geom_tiplab(d, aes(label=pond))

plot(njo, direction="upwards")

ggtree(njo) %<+% d +
geom_tiplab(aes(label=A)) +
geom_tiplab(aes(label=C)) +
geom_tiplab(aes(label=pond, color=pond), offset=1, angle=90, align=T) +
coord_flip()



#### demo
### phased
  seqSetFilter(genofile, variant.id=1:10, sample.id=c(
  "Spring_2017_DBunk_113",
  "Spring_2017_DBunk_122",
  "Spring_2017_DBunk_125",
  "Spring_2017_DBunk_143"))


  tmp <- as.data.table(getGenotype(genofile))
  tmp[,sample.id:=seqGetData(genofile, "sample.id")]

  dat.phase <- melt(tmp, id.vars="sample.id", variable.name="variant.id", value.name="geno")

  dat.phase[,allele1:=tstrsplit(geno, "\\|")[[1]]]
  dat.phase[,allele2:=tstrsplit(geno, "\\|")[[2]]]

### convert back to wide
  mat <- t(as.matrix(dcast(dat.phase, variant.id~sample.id, value.var=c("allele1", "allele2"))))
  dimnames(mat)[[2]] <- mat[1,]
  mat <- (mat[-1,])
  matn <- apply(mat, 2, as.numeric)

  dimnames(matn)[[1]] <- dimnames(mat)[[1]]
  dimnames(matn)[[2]] <- dimnames(mat)[[2]]

### nj tree
  dist.mat <- dist(matn)
  njo <- nj(dist.mat)
