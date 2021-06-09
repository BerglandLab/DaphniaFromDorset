#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(SeqArray)
  library(data.table)
  library(foreach)
  library(SeqVarTools)
  library(ape)

### open gds
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

### libraries
  library(data.table)
  library(ape)
  library(ggtree)
  library(ggplot2)

  load("~/mats.Rdata")

### iterate through
  foreach(i=1:length(mats))%do%{
    #i<-3
    matn <- mats[[i]]

    dm <- dist(matn, diag=T, upper=T)
    njo <- nj(dm)
    njo <- root(njo, outgroup="allele1_Spring_2016_W6_6.3")
    d <- data.table(label=njo$tip.label)
    d[,pond:=tstrsplit(label, "_")[[4]]]
    d[,A:=ifelse(grepl("April_2017_D8_213", label), "A", "")]
    d[,C:=ifelse(grepl("April_2017_D8_151", label), "C", "")]
    d <- as.data.frame(d)


    p <- ggtree(njo) %<+% d +
    geom_tiplab(aes(label=A)) +
    geom_tiplab(aes(label=C)) +
    geom_tiplab(aes(label=pond, color=pond), offset=1, angle=90, align=T) +
    coord_flip()

    ggsave(p, file=paste("~/qtl", i, ".pdf", sep=""))
  }



















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
