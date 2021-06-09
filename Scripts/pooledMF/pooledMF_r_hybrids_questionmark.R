#ijob -c1 -p standard -A berglandlab
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R


### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

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


  ac.inform <- ac.fd[(A.geno=="12" & C.geno=="11") |
                     (A.geno=="12" & C.geno=="22") |
                     (A.geno=="11" & C.geno=="12") |
                     (A.geno=="22" & C.geno=="12") |
                     (A.geno=="12" & C.geno=="12") |
                     (A.geno=="11" & C.geno=="22") |
                     (A.geno=="22" & C.geno=="11")]

  #ac.inform <- ac.inform[A.delta < 0.05 & C.delta < 0.05]
  #save(ac.inform, file="/scratch/aob2x/daphnia_hwe_sims/ac_inform.Rdata")


### pooled data
  load("/nv/vol186/bergland-lab/alan/totalADRDlongall.Rdata")
  geno <- geno[pond=="D8"]
  geno[,effRD:=floor(RRD)]
  geno[,effPA:=round(propalt*effRD)/effRD]

  geno.w <- dcast(geno[pond=="D8"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))
  geno.w[, f.hat := (effPA_D8Male1 + effPA_D8Male2 + effPA_D8PE1+ effPA_D8PE2)/4]
  geno.w[f.hat>0 & f.hat<1]


### save object
  save(geno, geno.w, ac.inform, snp.dt, file="~/pooled_f1_hybrids_questionmark.Rdata")

#### download
  scp aob2x@rivanna.hpc.virginia.edu:~/pooled_f1_hybrids_questionmark.Rdata ~/.
  R

### libraries
  library(data.table)
  library(ggplot2)

### load data
  load("~/pooled_f1_hybrids_questionmark.Rdata")

### basics
  hist(geno.w$f.hat, breaks=100)

  hist(geno.w[f.hat>0 & f.hat<1]$f.hat, breaks=1000)

### SNP that were retained in the individual based VCF file
  setkey(geno.w, chr, pos)
  setkey(snp.dt, chr, pos)
  geno.w.use <- merge(geno.w[f.hat>0 & f.hat<1], snp.dt)

### merge with "informative" A&C sites
  setkey(ac.inform, chr, pos)
  setkey(geno.w.use , chr, pos)

  dim(merge(ac.inform, geno.w.use))
  dim(geno.w.use)
  dim(ac.inform)

  m <- merge(ac.inform, geno.w.use, all.x=T, all.y=T)

### table
  prop.table(table(poolSeq= !is.na(m[f.hat>.05 & f.hat<.95]$effRD_D8Male1),
                   ACPoly = !is.na(m[f.hat>.05 & f.hat<.95]$A.geno))



  dim(ac.inform)
  dim(geno.w)
  dim(m)
### how many?
  table(poolSeq=!is.na(m$effRD_D8Male1), ACPoly=!is.na(m$A.geno))

### are pooled D8 samples basically F1s?
  m <- na.omit(m)
  m[,m.hat:=(effPA_D8Male1 * effRD_D8Male1 + effPA_D8Male2 * effRD_D8Male2) / (effRD_D8Male1 + effRD_D8Male2)]
  m[,f.hat:=(effPA_D8PE1 * effRD_D8PE1 + effPA_D8PE2 * effRD_D8PE2) / (effRD_D8PE1 + effRD_D8PE2)]
  m[!is.na(A),A.geno := unlist(sapply(m[!is.na(A)]$A, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[!is.na(B),B.geno := unlist(sapply(m[!is.na(B)]$B, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]



  m.ag <- m[,list(mu=c(mean(m.hat), mean(f.hat)),
          sd=c(sd(m.hat), sd(f.hat)),
          pool=c("m", "f")),
      list(A.geno, B.geno)]
  m.ag[,exp.fq:=(A.geno+B.geno)/2]


  save(m.ag, file="/nv/vol186/bergland-lab/alan/mf_expectation.Rdata")



  library(ggplot2)
  library(data.table)
  library(cowplot)

  load("/mnt/sammas_storage/bergland-lab/alan/mf_expectation.Rdata")

  F1.plot <- ggplot(data=m.ag, aes(x=pool, y=mu)) +
  geom_hline(aes(yintercept=exp.fq/2), linetype="solid", color="red", size=.75) +
  geom_errorbar(aes(ymin=mu-sd, ymax=mu+sd), width=.2, size=.75) +
  geom_point(size=4, fill="white", shape=21, stroke=1) +
  facet_grid(~A.geno+B.geno) +
  ylab("Frequencuy") +
  xlab("Pool") +
  scale_x_discrete(labels=c("m" = "Male", "f" = "Female")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave("~/F1_plot.pdf")
