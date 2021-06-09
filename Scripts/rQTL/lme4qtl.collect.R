## ijob -c1 -p largemem --mem 5G -A berglandlab

# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(foreach)
  library(data.table)
  library(doMC)
  registerDoMC(20)

### param
  #set <- "AxC"
  #set <- "CxC"
  set <- "all"

  ## /project/berglandlab/alan/lme4qtl/


### load pio.unique.n
  load(file=paste("/scratch/aob2x/daphnia_hwe_sims/lmer4qtl/pio.uniq.set.", set, ".Rdata", sep=""))
  pio.uniq.n.long <- pio.uniq.n[,list(setId=as.numeric(unlist(tstrsplit(ids, ";"))),
                                      setPos=as.numeric(unlist(tstrsplit(poss, ";"))),
                                      setChr=unlist(tstrsplit(chrs, ";"))),
                                list(id)]
  setkey(pio.uniq.n.long, id)



### load data
    fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/lmer4qtl/", set, full.name=T)
    fn <- fn[grepl("v3", fn)]
    o <- foreach(fn.i=fn)%dopar%{
      #fn.i <- fn[1]
      message(fn.i)
      load(fn.i)

      setkey(lmer.gwas, id)

      lmer.gwas.long <- merge(lmer.gwas, pio.uniq.n.long, allow.cartesian=T)

      setnames(lmer.gwas.long, c("chr", "pos", "id", "setId", "setPos", "setChr"),
                               c("marker.chr", "marker.pos", "marker.id", "id", "pos", "chr"))
      return(lmer.gwas.long)
    }
    o <- rbindlist(o)

  #o.ag <- o[, list(p.z=mean(p.z, na.rm=T), chisq=max(chisq), p.aov=min(p.aov, na.rm=T)), list(term, perm, id, chr, pos)]
  #o.ag.ag <- o.ag[,list(pr=mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id, chr, pos)]
  #o.ag.ag[pr==1, pr:=1/201]
  o.obs <- o[[1]]


    save(o, file=paste("~/lme4qtl_output.v3.", set, ".long.Rdata", sep=""))
    save(o.obs, file=paste("~/lme4qtl_output.v3.", set, ".obs.long.Rdata", sep=""))

    save(o, file=paste("/project/berglandlab/alan/lme4qtl_output.v3.", set, ".long.Rdata", sep=""))
























### scp aob2x@rivanna.hpc.virginia.edu:~/lme4qtl_output.AxC.long.Rdata ~/.
  library(data.table)
  library(ggplot2)

  load("~/lme4qtl_output.AxC.long.Rdata")

  o.ag <- o[, list(p.z=mean(p.z, na.rm=T), chisq=max(chisq), p.aov=min(p.aov, na.rm=T)), list(term, perm, id, chr, pos)]
  o.ag.ag <- o.ag[,list(pr=mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id, chr, pos)]
  o.ag.ag[pr==1, pr:=1/201]

### ANOVA output
  o.ag.perm <- o.ag[perm!=0, list(p.aov=c(quantile(p.aov, .01), quantile(p.aov, .05)),
                                  q=c(.01, .05)), list(term, chr)]

  ggplot() +
  geom_line(data=o.ag[perm==0], aes(x=pos, y=-log10(p.aov), color=chr)) +
  geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.aov), linetype=as.factor(q))) +
  facet_grid(term~chr)


### beta-p output
  o.ag.perm <- o.ag[perm!=0, list(p.z=quantile(p.z, .05)), list(term, chr)]

  ggplot() +
  geom_line(data=o.ag[perm==0], aes(x=pos, y=-log10(p.z), color=chr)) +
  geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.z))) +
  facet_grid(term~chr)

### normalized
  ggplot() +
  geom_line(data=o.ag.ag, aes(x=pos, y=-log10(1-pr), color=chr)) +
  facet_grid(term~chr)



### summarize




  o.ag.ag <- o.ag[,list(pr=mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id)]

  o.ag[perm>0, list(q=quantile(minp, .05)), list(term)]



### libraries
  library(SeqArray)

  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile))
  setkey(snp.dt, id)

  snp.dt[id==976200]
  seqSetFilter(genofile, variant.id=976200)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")


  target <- data.table(dosage=seqGetData(genofile, "$dosage"),
                      clone=seqGetData(genofile, "sample.id"))

  target <- merge(target, sc, "clone")

  table(na.omit(target[Species=="pulex"][population=="D8"][,list(dosage=mean(dosage.V1, na.rm=T)), list(SCnum)])$dosage)


)



ggplot() +
geom_line(data=gprime, aes(x=POS, y=G, color=CHROM)) +
geom_point(data=gprime[POS==3968943], aes(x=POS, y=G)) +
facet_grid(~CHROM)
