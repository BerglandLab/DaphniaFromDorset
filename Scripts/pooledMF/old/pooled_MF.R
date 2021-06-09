
### librarioes
  library(data.table)
  library(foreach)

### get job ID
  args <- commandArgs(trailingOnly = TRUE)

  nJobs <- as.numeric(args[1])
  jobId <- as.numeric(args[2])

  ### nJobs=1000; jobId=100
  print(nJobs)
  print(jobId)

### load data
  load("/project/berglandlab/Karen/SingleMomsMales20182019/PooledMomsMales/totrdfilt.Rdata")
  totrdfilt[,effRD:=floor(RRD)]
  totrdfilt[,effPA:=round(propalt*effRD)/effRD]


### make job array table
  snps.dt <- totrdfilt[,list(pos=unique(pos)), list(chr)]
  snps.dt[,id:=c(1:dim(snps.dt)[1])]
  groups <-  split(snps.dt$id, ceiling(seq_along(snps.dt$id)/(length(snps.dt$id)/nJobs)))

### get tmp data
  setkey(snps.dt, id)
  snps.dt.tmp <- snps.dt[J(groups[[jobId]])]

  setkey(totrdfilt, chr, pos)
  setkey(snps.dt.tmp, chr, pos)
  #totrdfilt.tmp <- totrdfilt[J(snps.dt.tmp)]

### define contrasts
pairs=list(c("D8Male1", "D8PE1"),
           c("D8Male1", "D8PE2"),
           c("D8Male2", "D8PE1"),
           c("D8Male2", "D8PE2"),
           c("D8Male1", "D8Male2"),
           c("D8PE1", "D8PE2"),
           c("DBunkMale", "DBunkPE1"),
           c("DBunkMale", "DBunkPE2"),
           c("DBunkPE1", "DBunkPE2"),
           c("DCatMale", "DCatPE1"),
           c("DCatMale", "DCatPE2"),
           c("DCatPE1", "DCatPE2"))


### iterate through SNPs



  o <- foreach(i=1:dim(snps.dt.tmp)[1], .combine="rbind", .errorhandling="remove")%do%{
    #i<-1
    print(paste(i, dim(snps.dt.tmp)[1], sep=" / "))
    tmp <- totrdfilt[J(snps.dt.tmp[i])]


    foreach(k=pairs, .combine="rbind", .errorhandling="remove")%do%{
      #k <- pairs[[1]]
      s1 <- tmp[Sample==k[1]]
      s2 <- tmp[Sample==k[2]]

      mat <- matrix(c(s1$effPA*s1$effRD, (1-s1$effPA)*s1$effRD,
                      s2$effPA*s2$effRD, (1-s2$effPA)*s2$effRD), nrow=2, byrow=T)

      fet <- fisher.test(mat)
      data.table(pairs=paste(k, collapse="_"),
                 chr=s1$chr, pos=s1$pos,
                 p=fet$p.value,
                 or=fet$estimate,
                 effPA.s1=s1$effPA, effRD.s1=s1$effRD,
                 effPA.s2=s2$effPA, effRD.s2=s2$effRD)
    }

   #t1 <- glm(effPA~sex + pond, data=tmp, weights=effRD, family=binomial())
   #t1.D8 <- glm(effPA~sex , data=tmp[pond=="D8"], weights=effRD, family=binomial())
   #t1.DBunk <- glm(effPA~sex , data=tmp[pond=="DBunk"], weights=effRD, family=binomial())
   #t1.DCat <- glm(effPA~sex , data=tmp[pond=="DCat"], weights=effRD, family=binomial())

   #data.table(chr=tmp$chr[1], pos=tmp$pos[1],
   #           p.value=c(summary(t1)$coef[2,4],
   #                     summary(t1.D8)$coef[2,4],
   #                     summary(t1.DBunk)$coef[2,4],
   #                     summary(t1.DCat)$coef[2,4]),
   #            mod=c("all", "D8", "DBunk", "DCat"))

  }

### save
  write.csv(o, file=paste("/scratch/aob2x/daphnia_hwe_sims/mf/mf_fet_", jobId, ".csv", sep=""), quote=F, row.names=F)
