# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  #library(qtl)
  library(data.table)
  library(lme4qtl)
  library(doMC)
  registerDoMC(10)

  #library(rrBLUP)
  #library(ggplot2)
  #library(patchwork)
  #library(lme4)


############################
#### Prep the input data ###
############################
  ### `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.convert_output.vitPath.R` does this

########################
### [P]henotype data ###
########################
  ### Load SC data
    sc <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020//Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  ## Male rates
    ### load raw male phenotype data
      males <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/MaleCensus.csv")
      #males <- fread("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rqtl.convertData.R")

    ### some renaming
     sconeliter <- sc[OneLiterPheno==1]
     sconeliter$SCB <- ifelse(sconeliter$LabGenerated==0 & sconeliter$SC!="OO", sconeliter$SC, sconeliter$clone)
     sconeliter$Clone <- sconeliter$clone

    ### tidy
     setkey(males, Clone)
     setkey(sconeliter, Clone)
     mmales <- merge(males, sconeliter)
     mmales$propmale <- mmales$Males/mmales$NewTotal
     mmales$propmalenoneo <- mmales$Males/(mmales$NewTotal-mmales$Neos)

     male <- mmales[,c("Clone", "Replicate", "Males", "NewTotal", "Neos", "propmale", "SC", "SCB", "AxCF1Hybrid", "OneLiterPheno"), with=F]
     setnames(male, "Clone", "clone")


  ## Epphipia fill rates and production
    ### load raw epphiphial fill data
     epp <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/EphippiaFinal.csv")
     sconeliter <- sc[OneLiterPheno==1]
     sconeliter$SCB <- ifelse(sconeliter$LabGenerated==0 & sconeliter$SC!="OO", sconeliter$SC, sconeliter$clone)
     sconeliter$Clone <- sconeliter$clone

     setkey(epp, Clone)
     setkey(sconeliter, Clone)
     mepp <- merge(epp, sconeliter)
     epp <- mepp

   # Add new variables
     epp$TotalEppB <- epp$TotalEpp*2
     epp$fill <- epp$TotalEmbCorr/epp$TotalEppB
     epp$fullid <- paste(epp$Clone,epp$Rep,sep="_")
     epp$Rep <- as.factor(epp$Rep)
     epp$Clone <- as.factor(epp$Clone)
     epp$Dayssince1stDate <- as.factor(epp$Dayssince1stDate)
     epp$SCB <- as.factor(epp$SCB)

     setkey(epp, "Clone")

### add in cross type
    epp[SCB=="A", gr:="A"]
    epp[SCB=="C", gr:="C"]
    epp[AxCF1Hybrid==1, gr:="AxC"]
    epp[OneLiterPheno==1 & AxCF1Hybrid==0 & SC=="selfedC", gr:="CxC"]

    table(epp$gr)

    male[SCB=="A", gr:="A"]
    male[SCB=="C", gr:="C"]
    male[AxCF1Hybrid==1, gr:="AxC"]
    male[OneLiterPheno==1 & AxCF1Hybrid==0 & SC=="selfedC", gr:="CxC"]

    table(epp$gr)

    table(male$gr)


#####################
### Genotype data ###
#####################

  ### load genotype data (from DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.convert_output.vitPath.R)
    pio <- fread(file= "/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all.vitPath.pio.csv")
    pio[,phase.fold.geno:=phase.geno]
    pio[phase.fold.geno%in%c(12, 21), phase.fold.geno:=12]


  ### identify unique SNPs (there shoudl be about 5K of them)
    pio.unique <- pio[,list(vec=paste(phase.fold.geno, collapse=".")), list(chr, pos, id)]
    pio.uniq.tab <- table(pio.unique$vec)
    pio.uniq.n <- pio.unique[,list(id=id[1], chr=chr[1], pos=pos[1], .N), list(vec)]

    setkey(pio, id, chr, pos)
    setkey(pio.uniq.n, id, chr, pos)
    pio.u <- merge(pio, pio.uniq.n[,-c("vec"), with=F])


    ((table(pio.u$imputedGeno, pio.u$phase.fold.geno)))
    ((table(pio.u$obs.dosage, pio.u$phase.fold.geno)))


  ### make PIO unique into wide
      pio.u.ag <-   pio.u[,list(phase.fold.geno=phase.fold.geno[1],
                            ng=length(unique(phase.fold.geno))),
                      list(chr, pos, id, clone)]

      pio.u.ag[,d:=as.numeric(factor(phase.fold.geno)) - 2]
      pio.u.ag.w <- dcast(pio.u.ag, clone ~ id, value.var="d")
      pio.u.ag.ag <- pio.u.ag[,list(chr=chr[1], pos=pos[1]), list(id)]


### make kinship matrix using the rrBLUP package

  kinship_matrix <- rrBLUP::A.mat(pio.u.ag.w[,-c("clone"),with=F], n.core = 10)
  dimnames(kinship_matrix)[[1]] <- (pio.u.ag.w$clone)
  dimnames(kinship_matrix)[[2]] <- (pio.u.ag.w$clone)

#### using lme4qtl
  perm <- 0
  lmer.gwas <- foreach(i=pio.ag.ag$id[3000:3050], .errorhandling="remove")%dopar%{

    message(which(i==pio.ag.ag$id))
    #i=2920958
    geno.tmp <- pio.u[id==i]

    ### male
      #perm <- 0
      if(perm==0) {
        gp.male.tmp <- merge(male, geno.tmp, by="clone")

      } else if(perm>0) {

        male.ag <- male[,list(.N), list(clone)]
        male.ag[,permClone:=sample(clone)]

        male.perm <- merge(male, male.ag, by="clone")
        setnames(male.perm, c("clone", "permClone"), c("origClone", "clone"))

        gp.male.tmp <- merge(male.perm, geno.tmp, by="clone")

      }

      gp.male.tmp[,d:=as.numeric(as.factor(phase.fold.geno))-2]

      male.full <- relmatGlmer(cbind(Males, NewTotal-Males-Neos) ~ as.factor(d) + (1|clone) + (1|Replicate:clone),
                          family=binomial(),
                          data=gp.male.tmp,
                          relmat=list(clone=kinship_matrix))

      male.red <- relmatGlmer(cbind(Males, NewTotal-Males-Neos) ~ (1|clone) + (1|Replicate:clone),
                          family=binomial(),
                          data=gp.male.tmp,
                          relmat=list(clone=kinship_matrix))

      aov.male <- anova(male.full, male.red)



    ### fill rate
      if(perm==0) {
        gp.fill.tmp <- merge(epp, geno.tmp, by="clone")

      } else if (perm>0) {

        epp.ag <- epp[,list(.N), list(clone)]
        epp.ag[,permClone:=sample(clone)]

        epp.perm <- merge(epp, epp.ag, by="clone")
        setnames(epp.perm, c("clone", "permClone"), c("origClone", "clone"))

        gp.fill.tmp <- merge(epp.perm, geno.tmp, by="clone")

      }

      gp.fill.tmp[,d:=as.numeric(as.factor(phase.fold.geno))-2]
      epp.full <- relmatGlmer(cbind(fill*TotalEppB, (1-fill)*TotalEppB) ~ as.factor(d) + (1|clone) + (1|Rep:clone),
                          family=binomial(),
                          data=gp.fill.tmp,
                          relmat=list(clone=kinship_matrix))

      epp.red <- relmatGlmer(cbind(fill*TotalEppB, (1-fill)*TotalEppB) ~ (1|clone) + (1|Rep:clone),
                          family=binomial(),
                          data=gp.fill.tmp,
                          relmat=list(clone=kinship_matrix))

      aov.fill <- anova(epp.full, epp.red)

    ### output
      averages <- merge(gp.male.tmp[,list(males.mu=sum(Males)/sum(NewTotal-Neos), males.n=sum(NewTotal-Neos)), list(phase.fold.geno)],
                        gp.fill.tmp[,list(fill.mu=sum(fill*TotalEppB, na.rm=T)/sum((1-fill)*TotalEppB, na.rm=T), fill.n=sum(TotalEppB, na.rm=T)), list(phase.fold.geno)],
                        by="phase.fold.geno")
      averages[,id:=i]

      out <- data.table(id=i,
                  chisq=c(aov.fill[2,6], aov.male[2,6]),
                  p=c(aov.fill[2,8], aov.male[2,8]),
                  term=c("fill", "male"))
      out[,perm:=perm]
      merge(out, averages, by="id", allow.cartesian=T)

  }
  lmer.gwas <- rbindlist(lmer.gwas)




























  ### BLUP
    male.glmer <- glmer(cbind(Males, NewTotal-Neos-Males) ~ gr + (1|Replicate)  + (1|clone), male, family=binomial())
    male.glmer2 <- glmer(cbind(Males, NewTotal-Neos-Males) ~ gr + (1|Replicate), male, family=binomial())
    anova(male.glmer, male.glmer2)

    epp.glmer <- glmer(cbind(fill*TotalEppB, (1-fill)*TotalEppB) ~ gr + (1|Rep)  + (1|clone), epp, family=binomial())
    epp.glmer2 <- glmer(cbind(fill*TotalEppB, (1-fill)*TotalEppB) ~ gr + (1|Rep), epp, family=binomial())
    anova(epp.glmer, epp.glmer2)

    blup <- merge(
                  data.table(male.blup=ranef(male.glmer)$clone[,1],
                            clone=rownames(ranef(male.glmer)$clone)),
                  data.table(epp.blup=ranef(epp.glmer)$clone[,1],
                            clone=rownames(ranef(epp.glmer)$clone)),
                  all=F)

  ### averages
    epp.ag <- epp[,list(emb=sum(fill*TotalEppB, na.rm=T),
                                  Nemb=sum(TotalEppB)),
                    list(SCB, clone, gr)]
    epp.ag[,tProp.emb:=asin(sqrt(emb/Nemb))]
    epp.ag <- na.omit(epp.ag)


    male.ag <- male[,list(Males=sum(Males),
                                  N=sum(NewTotal - Neos)),
                             list(SCB, clone, gr)]
    male.ag[,tProp.male:=asin(sqrt(Males/N))]
    male.ag <- na.omit(male.ag)

    pheno.ag <- merge(epp.ag, male.ag, by="clone")


library(ggplot2)
library(data.table)
load("~/foo.Rdata")


ggplot(data=lmer.gwas, aes(x=id, y=-log10(fill))) + geom_line()










# id.i <- 2142150 - max Gprime peak


library(bWGR)


epp.ag <- epp[,list(emb=sum(fill*TotalEppB, na.rm=T),
                              Nemb=sum(TotalEppB)),
                list(SCB, clone, gr)]
epp.ag[,tProp.emb:=asin(sqrt(emb/Nemb))]
epp.ag <- na.omit(epp.ag)


male.ag <- male[,list(Males=sum(Males),
                              N=sum(NewTotal - Neos)),
                         list(SCB, clone, gr)]
male.ag[,tProp.male:=asin(sqrt(Males/N))]
male.ag <- na.omit(male.ag)

pheno.ag <- merge(epp.ag, male.ag, by="clone")


pio.ag <-   pio[,list(phase.fold.geno=phase.fold.geno[1],
                      ng=length(unique(phase.fold.geno))),
                list(chr, pos, id, clone)]

pio.ag[,d:=as.numeric(factor(phase.fold.geno))]
pio.ag.w <- dcast(pio.ag, clone ~ id, value.var="d")

setkey(pio.ag.w, clone)
setkey(pheno.ag, clone)




pio.ag.w.pheno <- merge(pheno.ag, pio.ag.w)

names(pio.ag.w.pheno)[1:12]
y <- pio.ag.w.pheno$tProp

fit_BayesL = wgr(y=pio.ag.w.pheno$tProp,
                 X=as.matrix(pio.ag.w.pheno[,-c("clone","SCB","gr","Males","N","tProp"),with=F]),
                 de=TRUE)
cor(pio.ag.w.pheno$tProp,fit_BayesL$hat)


fit_BayesB <- BayesB(y=pio.ag.w.pheno$tProp,
                     X=as.matrix(pio.ag.w.pheno[,-c("clone","SCB.x","SCB.y""gr.x","gr.y", "emb", "Nemb", "tProp.emb", "Males", "N", "tProp.male"),with=F]))

cor(pio.ag.w.pheno$tProp, fit_BayesB$hat)




pio.ag.ag <- pio.ag[,list(vec=paste(d, collapse="")), list(chr, pos)]
setkey(pio.ag.ag, vec)
pio.ag.ag[,chrpos:=paste(chr, pos, sep=";")]
pio.unique <- pio.ag.ag[,list(chrpos=chrpos[1], n=length(chrpos)), list(vec, chr)]
pio.unique[,pos:=as.numeric(tstrsplit(chrpos, ";")[[2]])]

setkey(pio.unique, chr, pos)
setkey(pio.ag, chr, pos)

pio.ag.w <- dcast(pio.ag[J(pio.unique[,c("chr", "pos"),with=F])], clone ~ id, value.var="d")

setkey(pio.ag.w, clone)
setkey(pheno.ag, clone)




pio.ag.w.pheno <- merge(pheno.ag, pio.ag.w)







#### statgetGWAS

  pio.ag.ag <- pio.ag[,list(vec=paste(d, collapse="")), list(chr, pos)]
  setkey(pio.ag.ag, vec)
  pio.ag.ag[,chrpos:=paste(chr, pos, sep=";")]
  pio.unique <- pio.ag.ag[,list(chrpos=chrpos[1], n=length(chrpos)), list(vec, chr)]
  pio.unique[,pos:=as.numeric(tstrsplit(chrpos, ";")[[2]])]

  setkey(pio.unique, chr, pos)
  setkey(pio.ag, chr, pos)

  pio.ag.w <- dcast(pio.ag[J(pio.unique[,c("chr", "pos"),with=F])], clone ~ id, value.var="d")

  setkey(pio.ag.w, clone)
  setkey(pheno.ag, clone)




  pio.ag.w.pheno <- merge(pheno.ag, pio.ag.w)



X.bm <- as.matrix(pio.ag.w.pheno[,-c("clone","male.blup", "epp.blup"),with=F])-1
row.names(X.bm) <- pio.ag.w.pheno$clone
dimnames(X.bm)[[2]] <- names(pio.ag.w.pheno[,-c("clone","male.blup", "epp.blup"),with=F])

map <- as.data.frame(pio.ag.ag[,c("chr", "pos"), with=F])
rownames(map) <- pio.ag.ag$id
gDataDrops <- createGData(geno = X.bm, map = map)

pheno <- data.table(pio.ag.w.pheno[,c("clone", "male.blup", "epp.blup"),with=F])
names(pheno)[1] <- "genotype"
pheno$genotype <- as.factor(pheno$genotype)

gDataDrops <- createGData(gData = gDataDrops, pheno = as.data.frame(pheno))


summary(gDataDrops)
gDataDropsDedup <- codeMarkers(gDataDrops, impute = FALSE, verbose = TRUE, removeDuplicates = F)

GWASDrops <- runSingleTraitGwas(gData = gDataDropsDedup,
                                traits = c("male.blup", "epp.blup"),
                                kinshipMethod="vanRaden", genomicControl=F)

summary(GWASDrops$GWAResult[[1]])
summary(GWASDrops)

tmp <- as.data.table(GWASDrops$GWAResult)
tmp[as.data.frame.pheno..LOD>4]


library(biglasso)

### lasso
X.bm <- as.big.matrix(as.matrix(pio.ag.w.pheno[,-c("clone","SCB.x","SCB.y","gr.x","gr.y", "emb", "Nemb", "tProp.emb", "Males", "N", "tProp.male"),with=F]))
X.bm <- as.big.matrix(as.matrix(pio.ag.w.pheno[,-c("clone","male.blup", "epp.blup"),with=F]))

### run lasso models
if(!file.exists(pred.out.fn)) {
message("running lasso")


### embryo fill rate
  cvfit.emb <- cv.biglasso(X=X.bm, y=pio.ag.w.pheno$tProp.emb,
                          nfolds = 10, ncores = 10,
                          family="gaussian",
                          screen="SSR-Slores",
                          penalty="enet",
                          max.iter=500000,
                          alpha=.5,
                          trace=T)
  summary(cvfit.emb)


  coefs.tmp <- coef(cvfit.emb)

  coef.tmp <- data.table(coef=coefs.tmp[which(coefs.tmp != 0)],
                        id=as.numeric(dimnames(coefs.tmp)[[1]][which(coefs.tmp != 0)]))

  pio.ag.ag <- pio.ag[,list(.N), list(chr, pos, id)]

  lasso.dt <- merge(coef.tmp, pio.ag.ag, by="id")
  lasso.dt[,.N, chr]

### male production
  cvfit.male <- cv.biglasso(X=X.bm, y=pio.ag.w.pheno$tProp.male,
                          ncores = 10,
                          nfolds = 10,
                          family="gaussian",
                          screen="SSR-Slores",
                          penalty="enet",
                          max.iter=50000,
                          alpha=.5)
  summary(cvfit.male)


  coefs.tmp <- coef(cvfit.male)

  coef.tmp <- data.table(coef=coefs.tmp[which(coefs.tmp != 0)],
                        id=as.numeric(dimnames(coefs.tmp)[[1]][which(coefs.tmp != 0)]))

  pio.ag.ag <- pio.ag[,list(.N), list(chr, pos, id)]

  lasso.dt <- merge(coef.tmp, pio.ag.ag, by="id")
  lasso.dt[,.N, chr]


load("~/peaks.Rdata")

gprime <- as.data.table(gprime)
setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))

#gp <- ggplot(data=gprime, aes(x=pos, y=Gprime, color=chr)) + geom_line() + facet_grid(.~chr)
#lasso <- ggplot(data=lasso.dt, aes(x=pos, y=coef, color=chr)) + geom_point() + facet_grid(.~chr)

setkey(lasso.dt, chr, pos)
setkey(gprime, chr, pos)
m <- merge(lasso.dt, gprime, all=T)
m[Gprime>5][!is.na(coef)]


fisher.test(table(m$Gprime>15, !is.na(m$coef)))

save(m, file="~/peaks_lasso.Rdata")


scp aob2x@rivanna.hpc.virginia.edu:~/peaks_lasso.Rdata .






save(lasso.dt, file="~/lassoDT.Rdata")



scp aob2x@rivanna.hpc.virginia.edu:~/lassoDT.Rdata ~/.

library(data.table)
library(ggplot2)

load("~/peaks_lasso.Rdata")

ggplot() +
geom_line(data=m, aes(x=pos, y=Gprime, color=chr)) +
geom_point(data=m[!is.na(coef)], aes(x=pos, y=Gprime, color=chr)) +
facet_grid(.~chr)




lasso <- ggplot(data=lasso.dt, aes(x=pos, y=coef, color=chr)) + geom_point() + facet_grid(.~chr)

setkey(lasso.dt, chr, pos)
setkey(gprime, chr, pos)
m <- merge(lasso.dt, gprime, all=T)



library(patchwork)
gp/lasso







bl <- data.table(beta=fit_BayesL$b,
                id=names(pio.ag.w.pheno[,-c("clone","SCB","gr","emb","N","tProp"),with=F]))

gwa <- emGWA(y=pio.ag.w.pheno$tProp,
            gen=as.matrix(pio.ag.w.pheno[,-c("clone","SCB","gr","emb","N","tProp"),with=F]))
gwa.dt <- data.table(beta=gwa$b, p=gwa$PVAL,
                     id=names(pio.ag.w.pheno[,-c("clone","SCB","gr","emb","N","tProp"),with=F]))


h2 <- emCV(y=pio.ag.w.pheno$tProp,
            gen=as.matrix(pio.ag.w.pheno[,-c("clone","SCB","gr","emb","N","tProp"),with=F]))


pio.ag.ag <- pio.ag[,list(.N), list(chr, pos, id)]
gwa.dt[,id:=as.numeric(id)]

gwa.dt <- merge(gwa.dt, pio.ag.ag, by="id")

save(gwa.dt, file="~/gwaDT.Rdata")



### scp aob2x@rivanna.hpc.virginia.edu:~/bl.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:~/gwaDT.Rdata ~/.

library(ggplot2)
library(data.table)

load("~/bl.Rdata")
plot(beta~id, bl)

load("~/gwaDT.Rdata")


ggplot(data=gwa.dt, aes(x=pos, y=-log10(p), color=chr)) +
geom_point() +
facet_grid(~chr)


tmp.fill.ag <- tmp.fill[,list(emb=sum(fill*TotalEppB, na.rm=T),
                              N=sum(TotalEppB),
                              phase.fold.geno=phase.fold.geno[1],
                              ng=length(unique(phase.fold.geno))),
                        list(SCB, gr)]
tmp.fill.ag[,perm:=sample(phase.fold.geno)]


id.i <- pio.ag[chr=="Scaffold_9199_HRSCAF_10755" & pos==6232363]$id
pio.ag[,i:=c(1:dim(pio.ag)[1])]
pio.ag[i>=(70654-50) & i<=(70654+50)]


tmp.male <- merge(male, pio[id>=(id.i-100) & id<=(id.i+100)], by="clone")























scp aob2x@rivanna.hpc.virginia.edu:~/gwasOutput.Rdata ~/.
library(patchwork)
library(data.table)
library(ggplot2)
load("~/gwasOutput.Rdata")
load("~/peaks.Rdata")

setkey(out, chr, pos)
setkey(gprime, chr, pos)
m <- merge(out, gprime)
fisher.test(table(m$pvalue<.005, m$p<.005, m$term)[,,2])



gprime <- as.data.table(gprime)

id.i <- sapply(c(1:12), function(x) out[chr==peaks$CHROM[x]][which.min(abs(peaks$posMaxGprime[x] - pos))]$id)

id.i <- c(out[chr==peaks$CHROM[11]][which.min(abs(peaks$posMaxGprime[11] - pos))]$id,
          out[chr==peaks$CHROM[11]][which.min(abs(peaks$start[11] - pos))]$id,
          out[chr==peaks$CHROM[11]][which.min(abs(peaks$end[11] - pos))]$id)

out[id==id.i[6]]


qtl <- ggplot() +
geom_point(data=out, aes(x=pos, y=-log10(p), color=chr)) +
facet_grid(term~chr, scales="free") +
theme(strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +



#setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))
Gprime.plot <- ggplot(data=gprime, aes(x=pos, y=Gprime, color=chr)) +
#geom_rect(data=peaks, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="grey", alpha=.5) +
geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
geom_line(size=.75) +
facet_grid(.~chr, scales="free_x") +
theme(strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
ggtitle("pooledWild")

qtl / Gprime.plot








ggplot(data=out, aes(x=id, y=-log10(p))) + geom_point() + facet_wrap(~term)

### test
load(file="/nv/vol186/bergland-lab/alan/peaks.Rdata")

pio[chr==peaks$CHROM[11]][which.min(abs(peaks$posMaxGprime[11] - pos))]
tmp.male <- merge(mmales, pio[id==3042479], by="clone")


t1.male <- glmer(propmalenoneo ~ as.factor(phase.geno) + (1|SCB), tmp.male, family=binomial(),
                  weights=tmp.male$NewTotal - tmp.male$Neos)
summary(t1.male)


anova(lm(propmalenoneo~as.factor(phase.geno), tmp))
t1.male <- lm(propmalenoneo~as.factor(imputedGeno), tmp)





### simple gwas
pio.ag <- pio[,list(.N), list(chr, pos, id)]

ii <- seq(from=1, to=dim(pio.ag)[1], by=1)


  out <- foreach(i=ii, .errorhandling="remove")%dopar%{
    message(paste(which(i==ii), length(ii), sep=" / "))

    # id.i <- 2142150 - max Gprime peak
    # i <- 5158
    # id.i <- 3042561
    id.i <- pio.ag$id[i]

    ### male model: phase.geno
      tmp.male <- merge(male, pio[id==id.i], by="clone")
      tmp.male.ag <- tmp.male[,list(Males=sum(Males),
                                    N=sum(NewTotal - Neos), phase.fold.geno=phase.fold.geno[1], ng=length(unique(phase.fold.geno))),
                               list(SCB, gr)]

      #t1.male <- glmer(propmalenoneo ~ as.factor(phase.fold.geno) + (1|SCB) + (1|Replicate:SCB), tmp.male, family=binomial(),
      #                weights=tmp.male$NewTotal - tmp.male$Neos)
      #t0.male <- glmer(propmalenoneo ~ 1 + (1|SCB) + (1|Replicate:SCB), tmp.male, family=binomial(),
      #                weights=tmp.male$NewTotal - tmp.male$Neos)

    t1.male <- glm(cbind(Males, N-Males) ~ as.factor(phase.fold.geno) + gr, tmp.male.ag, family=binomial())
    t0.male <- glm(cbind(Males, N-Males) ~ gr, tmp.male.ag , family=binomial())

    aov.male <- anova(t1.male, t0.male, test="Chisq")

  ### fill rate model
    tmp.fill <- merge(epp, pio[id==id.i], by="clone")
    tmp.fill.ag <- tmp.fill[,list(emb=sum(fill*TotalEppB, na.rm=T),
                                  N=sum(TotalEppB),
                                  phase.fold.geno=phase.fold.geno[1],
                                  ng=length(unique(phase.fold.geno))),
                            list(SCB, gr)]
    tmp.fill.ag[,perm:=sample(phase.fold.geno)]



    #t1.fill <- glmer(fill ~ as.factor(phase.fold.geno) + (1|SCB) + (1|Rep:SCB), tmp.fill, family=binomial(),
    #                weights=tmp.fill$TotalEppB)
    #t0.fill <- glmer(fill ~ 1 + (1|SCB) + (1|SCB:Rep), tmp.fill, family=binomial(),
    #                weights=tmp.fill$TotalEppB)
    #anova(t1.fill, t0.fill, test="Chisq")

    t1.fill <- glm(cbind(emb, N-emb) ~ as.factor(phase.fold.geno) + gr , tmp.fill.ag, family=binomial())
    t0.fill <- glm(cbind(emb, N-emb) ~ gr , tmp.fill.ag, family=binomial())
    aov.fill <- anova(t1.fill, t0.fill, test="Chisq")

    #tmp.fill.ag[,perm:=sample(phase.fold.geno)]
    #t1.fill <- glm(cbind(emb, N-emb) ~ as.factor(perm) +gr , tmp.fill.ag, family=binomial())
    #t0.fill <- glm(cbind(emb, N-emb) ~ gr , tmp.fill.ag, family=binomial())
    #anova(t1.fill, t0.fill, test="Chisq")

    averages <- merge(tmp.male.ag[,list(males.mu=sum(Males)/sum(N), males.n=sum(N), geno.n=length(N)), list(phase.fold.geno)],
                      tmp.fill.ag[,list(fill.mu=sum(emb)/sum(N), fill.n=sum(N), geno.n=length(N)), list(phase.fold.geno)],
                      by="phase.fold.geno")
    averages[,id:=id.i]

    out <- data.table(id=id.i,
                chisq=c(aov.fill[2,4], aov.male[2,4]),
                p=c(aov.fill[2,5], aov.male[2,5]),
                term=c("fill", "male"))
    merge(out, averages, by="id", allow.cartesian=T)

  }

out <- rbindlist(out)
out[,pa:=p.adjust(p, "fdr")]
setkey(out, id)
setkey(pio.ag, id)

out <- merge(out, pio.ag)
save(out, file="~/gwasOutput.Rdata")
out[term=="fill"][which.min(chisq)]
out[id==977045]

out[id==3065515]
