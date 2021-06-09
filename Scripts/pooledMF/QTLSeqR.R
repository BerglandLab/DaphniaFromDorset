### ijob -c1 -p standard -A berglandlab
### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


#install.packages("devtools")
#
## use devtools to install QTLseqr
#devtools::install_github("bmansfeld/QTLseqr")


### libraries
  library(data.table)
  library(QTLseqr)
  library(foreach)

### load AC inform object made by `informative_sites.R`
  load(file="/scratch/aob2x/daphnia_hwe_sims/ac_inform.Rdata")
  ac.inform[!(A.geno==12 & C.geno)]

### 1. using the ASE read counter data and running gprime separately for each rep PE1 vs Male1; PE2 vs Male2
### load data
  gprime.peaks <- foreach(r=list(c(1,1), c(2,2), c(1,2), c(2,1)))%do%{
    message(paste(r, collapse="/"))
    male <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8Male.pooledAF.aseReadCounter.allvariant.Male", r[1], ".delim", sep=""))
    pe <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8PE.pooledAF.aseReadCounter.allvariant.PE", r[2], ".delim", sep=""))

    setnames(male,
             c("contig", "position", "refAllele", "altAllele", "refCount", "altCount"),
             c("CHROM", "POS", "REF", "ALT", "AD_REF.LOW", "AD_ALT.LOW"))

    male[,freq:=AD_ALT.LOW/(AD_ALT.LOW + AD_REF.LOW)]
    male[,rd:=AD_ALT.LOW + AD_REF.LOW]
    male[,effrd:=round((70*rd)/(70+rd))]

    male[,AD_REF.LOW:=round((1-freq)*effrd)]
    male[,AD_ALT.LOW:=round((freq)*effrd)]

   setnames(pe,
            c("contig", "position", "refAllele", "altAllele", "refCount", "altCount"),
            c("CHROM", "POS", "REF", "ALT", "AD_REF.HIGH", "AD_ALT.HIGH"))
    pe[,freq:=AD_ALT.HIGH/(AD_ALT.HIGH + AD_REF.HIGH)]
    pe[,rd:=AD_ALT.HIGH + AD_REF.HIGH]
    pe[,effrd:=round((100*rd)/(100+rd))]

    pe[,AD_REF.HIGH:=round((1-freq)*effrd)]
    pe[,AD_ALT.HIGH:=round((freq)*effrd)]


    setkey(male, CHROM, POS)
    setkey(pe, CHROM, POS)

    m <- merge(male[,c("CHROM", "POS", "REF", "ALT", "AD_REF.LOW", "AD_ALT.LOW"),with=F],
               pe[,c("CHROM", "POS", "REF", "ALT", "AD_REF.HIGH", "AD_ALT.HIGH"),with=F])

    m[,deltaSNP:=qlogis(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - qlogis(AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
    m.ag <- m[,.N,CHROM]
    m <-  m[CHROM%in%m.ag[N>1000]$CHROM]

    m[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
    m[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
    m[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]

    m <- m[!is.na(deltaSNP) & deltaSNP!=Inf & deltaSNP!=-Inf]
    ### filter?
    df_filt <- filterSNPs(SNPset = as.data.frame(m),
                         refAlleleFreq = 0.15,
                         minTotalDepth = 20,
                         maxTotalDepth = 5000,
                         depthDifference = 1500,
                         minSampleDepth = 20,
                         verbose = TRUE)


    gprime <- runGprimeAnalysis(SNPset=df_filt, windowSize=250000) #250000
    peaks <- getQTLTable(gprime, alpha=.05)
    peaks <- as.data.table(peaks)
    peaks[,rep:=paste(r, collapse="/")]

    gprime <- as.data.table(gprime)
    gprime[,rep:=paste(r, collapse="/")]

    m[,rep:=paste(r, collapse="/")]

    list(gprime, peaks, m)
  }
  gprime <- rbindlist(lapply(gprime.peaks, function(x) x[[1]]))
  peaks <- rbindlist(lapply(gprime.peaks, function(x) x[[2]]))
  alleleFreqs <- rbindlist(lapply(gprime.peaks, function(x) x[[3]]))

  table(peaks$rep)
  table(gprime$rep)

  save(gprime, peaks, file="/scratch/aob2x/daphnia_hwe_sims/gprime_peaks.replicates_combos.250K.05.Rdata")
  save(gprime, peaks, file="/project/berglandlab/alan/gprime_peaks.replicates_combos.250K.05.Rdata")
  save(alleleFreqs, file="/project/berglandlab/alan/alleleFreqs.replicates_combos.250K.05.Rdata")

  load(file="/project/berglandlab/alan/gprime_peaks.replicates.250K.05.Rdata")

### 2. using the ASE read counter data and running gprime separately for each rep PE1 vs PE2; Male1 vs Male2
### load data
  gprime.peaks.groups <- foreach(group=c("Male", "PE"))%do%{
    r1 <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8", group, ".pooledAF.aseReadCounter.allvariant.", group, "1.delim", sep=""))
    r2 <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8", group, ".pooledAF.aseReadCounter.allvariant.", group, "2.delim", sep=""))

    setnames(r1,
             c("contig", "position", "refAllele", "altAllele", "refCount", "altCount"),
             c("CHROM", "POS", "REF", "ALT", "AD_REF.LOW", "AD_ALT.LOW"))

    r1[,freq:=AD_ALT.LOW/(AD_ALT.LOW + AD_REF.LOW)]
    r1[,rd:=AD_ALT.LOW + AD_REF.LOW]
    r1[,effrd:=round((70*rd)/(70+rd))]

    r1[,AD_REF.LOW:=round((1-freq)*effrd)]
    r1[,AD_ALT.LOW:=round((freq)*effrd)]

   setnames(r2,
            c("contig", "position", "refAllele", "altAllele", "refCount", "altCount"),
            c("CHROM", "POS", "REF", "ALT", "AD_REF.HIGH", "AD_ALT.HIGH"))
    r2[,freq:=AD_ALT.HIGH/(AD_ALT.HIGH + AD_REF.HIGH)]
    r2[,rd:=AD_ALT.HIGH + AD_REF.HIGH]
    r2[,effrd:=round((100*rd)/(100+rd))]

    r2[,AD_REF.HIGH:=round((1-freq)*effrd)]
    r2[,AD_ALT.HIGH:=round((freq)*effrd)]


    setkey(r1, CHROM, POS)
    setkey(r2, CHROM, POS)

    m <- merge(r1[,c("CHROM", "POS", "REF", "ALT", "AD_REF.LOW", "AD_ALT.LOW"),with=F],
               r2[,c("CHROM", "POS", "REF", "ALT", "AD_REF.HIGH", "AD_ALT.HIGH"),with=F])

    m[,deltaSNP:=qlogis(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - qlogis(AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
    m.ag <- m[,.N,CHROM]
    m <-  m[CHROM%in%m.ag[N>1000]$CHROM]

    m[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
    m[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
    m[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]

    m <- m[!is.na(deltaSNP) & deltaSNP!=Inf & deltaSNP!=-Inf]
    ### filter?
    df_filt <- filterSNPs(SNPset = as.data.frame(m),
                         refAlleleFreq = 0.15,
                         minTotalDepth = 20,
                         maxTotalDepth = 5000,
                         depthDifference = 1500,
                         minSampleDepth = 20,
                         verbose = TRUE)


    gprime <- runGprimeAnalysis(SNPset=df_filt, windowSize=250000) #250000
    peaks <- getQTLTable(gprime, alpha=.05)
    peaks <- as.data.table(peaks)
    peaks[,group:=group]

    gprime <- as.data.table(gprime)
    gprime[,group:=group]

    m[,group:=group]

    list(gprime, peaks, m)
  }
  gprime.groups <- rbindlist(lapply(gprime.peaks.groups, function(x) x[[1]]))
  peaks.groups <- rbindlist(lapply(gprime.peaks.groups, function(x) x[[2]]))
  alleleFreqs.groups <- rbindlist(lapply(gprime.peaks.groups, function(x) x[[3]]))

  table(peaks$rep)

  save(gprime.groups, peaks.groups, file="/scratch/aob2x/daphnia_hwe_sims/gprime_peaks.groups.250K.05.Rdata")
  save(gprime.groups, peaks.groups, file="/project/berglandlab/alan/gprime_peaks.groups.250K.05.Rdata")
  save(alleleFreqs.groups, file="/project/berglandlab/alan/alleleFreqs.groups.250K.05.Rdata")


### 3. using the ASE read counter data and adding neff for ombibus test
  ### load data
    inputData <- foreach(r=c(1,2))%do%{
      #r<-1
      message(r)
      male <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8Male.pooledAF.aseReadCounter.allvariant.Male", r, ".delim", sep=""))
      pe <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8PE.pooledAF.aseReadCounter.allvariant.PE", r, ".delim", sep=""))

      setnames(male,
               c("contig", "position", "refAllele", "altAllele", "refCount",   "altCount"),
               c("CHROM",  "POS",      "REF",       "ALT",       "AD_REF.LOW", "AD_ALT.LOW"))

      male[,freq:=AD_ALT.LOW/(AD_ALT.LOW + AD_REF.LOW)]
      male[,rd:=AD_ALT.LOW + AD_REF.LOW]
      male[,effrd:=round((70*rd)/(70+rd))]

      male[,AD_REF.LOW.orig:=AD_REF.LOW]
      male[,AD_ALT.LOW.orig:=AD_ALT.LOW]

      male[,AD_REF.LOW.neff:=round((1-freq)*effrd)]
      male[,AD_ALT.LOW.neff:=round((freq)*effrd)]

     setnames(pe,
              c("contig", "position", "refAllele", "altAllele", "refCount",    "altCount"),
              c("CHROM",  "POS",      "REF",       "ALT",       "AD_REF.HIGH", "AD_ALT.HIGH"))
      pe[,freq:=AD_ALT.HIGH/(AD_ALT.HIGH + AD_REF.HIGH)]
      pe[,rd:=AD_ALT.HIGH + AD_REF.HIGH]
      pe[,effrd:=round((100*rd)/(100+rd))]

      pe[,AD_REF.HIGH.orig:=AD_REF.HIGH]
      pe[,AD_ALT.HIGH.orig:=AD_ALT.HIGH]

      pe[,AD_REF.HIGH.neff:=round((1-freq)*effrd)]
      pe[,AD_ALT.HIGH.neff:=round((freq)*effrd)]


      setkey(male, CHROM, POS)
      setkey(pe, CHROM, POS)

      m <- merge(male[,c("CHROM", "POS", "REF", "ALT", "AD_REF.LOW.orig", "AD_ALT.LOW.orig", "AD_REF.LOW.neff", "AD_ALT.LOW.neff", "effrd"),with=F],
                 pe[,c("CHROM", "POS", "REF", "ALT", "AD_REF.HIGH.orig", "AD_ALT.HIGH.orig", "AD_REF.HIGH.neff", "AD_ALT.HIGH.neff", "effrd"),with=F])

      m[,deltaSNP.orig:=qlogis(AD_ALT.LOW.orig/(AD_ALT.LOW.orig+AD_REF.LOW.orig)) - qlogis(AD_ALT.HIGH.orig/(AD_ALT.HIGH.orig+AD_REF.HIGH.orig))]
      m[,deltaSNP.neff:=qlogis(AD_ALT.LOW.neff/(AD_ALT.LOW.neff+AD_REF.LOW.neff)) - qlogis(AD_ALT.HIGH.neff/(AD_ALT.HIGH.neff+AD_REF.HIGH.neff))]

      m.ag <- m[,.N,CHROM]
      m <-  m[CHROM%in%m.ag[N>1000]$CHROM]

      m[,REF_FRQ.orig:=(AD_REF.LOW.orig + AD_REF.HIGH.orig) / (AD_REF.LOW.orig + AD_REF.HIGH.orig + AD_ALT.LOW.orig + AD_ALT.HIGH.orig )]
      m[,DP.LOW.orig:=  AD_REF.LOW.orig +  AD_ALT.LOW.orig]
      m[,DP.HIGH.orig:=AD_REF.HIGH.orig + AD_ALT.HIGH.orig]

      m[,REF_FRQ.neff:=(AD_REF.LOW.neff + AD_REF.HIGH.neff) / (AD_REF.LOW.neff + AD_REF.HIGH.neff + AD_ALT.LOW.neff + AD_ALT.HIGH.neff )]
      m[, DP.LOW.neff:= AD_REF.LOW.neff +  AD_ALT.LOW.neff]
      m[,DP.HIGH.neff:=AD_REF.HIGH.neff + AD_ALT.HIGH.neff]

      m[,rep:=r]
    }
    inputData <- rbindlist(inputData)

    save(inputData, file="/project/berglandlab/alan/PoolSeq_inputData.combined.Rdata")

  ### locally
   scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/PoolSeq_inputData.combined.Rdata ~/.

   load("~/PoolSeq_inputData.combined.Rdata")

   library(ggplot2)
   library(data.table)
   library(QTLseqr)


   m <- inputData[,
                   list(AD_REF.LOW=  sum(AD_REF.LOW.neff),
                        AD_ALT.LOW=  sum(AD_ALT.LOW.neff),
                        AD_REF.HIGH=sum(AD_REF.HIGH.neff),
                        AD_ALT.HIGH=sum(AD_ALT.HIGH.neff),
                        deltaSNP=     mean(deltaSNP.neff),
                        REF_FRQ=       mean(REF_FRQ.neff),
                        DP.LOW=          sum(DP.LOW.neff),
                        DP.HIGH=        sum(DP.HIGH.neff)),
                   list(CHROM, POS, REF=REF.x, ALT=ALT.x)]

   m <- m[deltaSNP!=-Inf & deltaSNP!=Inf]

   df_filt <- filterSNPs(SNPset = as.data.frame(m),
                        refAlleleFreq = 0.15,
                        minTotalDepth = 100,
                        maxTotalDepth = 5000,
                        depthDifference = 1500,
                        minSampleDepth = 10,
                        verbose = TRUE)

   gprime <- runGprimeAnalysis(SNPset=df_filt, windowSize=250000) #250000
   gprime <- as.data.table(gprime)
   #gprime[Gprime>40]

   peaks <- getQTLTable(gprime, alpha=.1)
   peaks <- as.data.table(peaks)
   peaks


   qtlseq <- runQTLseqAnalysis(SNPset=df_filt, windowSize=250000, bulkSize=40)
   #getQTLTable(SNPset=qtlseq, method="QTLseq")
   #negLog10Pval
   p <- plotQTLStats(SNPset = gprime, var="Gprime") +
        geom_hline(aes(yintercept=min(gprime[qvalue<.1]$Gprime)), color="red") +
        geom_hline(aes(yintercept=min(gprime[qvalue<.05]$Gprime)), color="blue") +
        geom_hline(aes(yintercept=min(gprime[qvalue<.01]$Gprime)), color="green")
  p


  load("/Users/alanbergland/lme4qtl_output.v3.all.obs.long.Rdata")


### scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/gprime_peaks.replicates.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/gprime_peaks.groups.250K.05.Rdata ~/.

  library(ggplot2)
  library(data.table)
  library(cowplot)
  load("~/gprime_peaks.replicates.Rdata")

  poolggplot() +
  geom_vline(data=peaks, aes(xintercept=posMaxGprime), color="grey") +
  geom_line(data=gprime, aes(x=POS, y=Gprime, color=CHROM)) +
  facet_grid(rep~CHROM) +
  theme(legend.position = "none")














































#### using the ASE read counter data and adding neff for ombibus test
  ### load data
    inputData <- foreach(r=c(1,2))%do%{
      #r<-1
      message(r)
      male <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8Male.pooledAF.aseReadCounter.allvariant.Male", r, ".delim", sep=""))
      pe <- fread(paste("/scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8PE.pooledAF.aseReadCounter.allvariant.PE", r, ".delim", sep=""))

      setnames(male,
               c("contig", "position", "refAllele", "altAllele", "refCount", "altCount"),
               c("CHROM", "POS", "REF", "ALT", "AD_REF.LOW", "AD_ALT.LOW"))

      male[,freq:=AD_ALT.LOW/(AD_ALT.LOW + AD_REF.LOW)]
      male[,rd:=AD_ALT.LOW + AD_REF.LOW]
      male[,effrd:=round((70*rd)/(70+rd))]

      male[,AD_REF.LOW:=round((1-freq)*effrd)]
      male[,AD_ALT.LOW:=round((freq)*effrd)]

     setnames(pe,
              c("contig", "position", "refAllele", "altAllele", "refCount", "altCount"),
              c("CHROM", "POS", "REF", "ALT", "AD_REF.HIGH", "AD_ALT.HIGH"))
      pe[,freq:=AD_ALT.HIGH/(AD_ALT.HIGH + AD_REF.HIGH)]
      pe[,rd:=AD_ALT.HIGH + AD_REF.HIGH]
      pe[,effrd:=round((100*rd)/(100+rd))]

      pe[,AD_REF.HIGH:=round((1-freq)*effrd)]
      pe[,AD_ALT.HIGH:=round((freq)*effrd)]


      setkey(male, CHROM, POS)
      setkey(pe, CHROM, POS)

      m <- merge(male[,c("CHROM", "POS", "REF", "ALT", "AD_REF.LOW", "AD_ALT.LOW", "effrd"),with=F],
                 pe[,c("CHROM", "POS", "REF", "ALT", "AD_REF.HIGH", "AD_ALT.HIGH", "effrd"),with=F])

      m[,deltaSNP:=qlogis(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - qlogis(AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
      m.ag <- m[,.N,CHROM]
      m <-  m[CHROM%in%m.ag[N>1000]$CHROM]

      m[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
      m[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
      m[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]
      m[,rep:=r]
    }
    inputData <- rbindlist(inputData)
    m <- inputData[,list(AD_REF.LOW=sum(AD_REF.LOW),
                         AD_ALT.LOW=sum(AD_ALT.LOW),
                         AD_REF.HIGH=sum(AD_REF.HIGH),
                         AD_ALT.HIGH=sum(AD_ALT.HIGH),
                         deltaSNP=mean(deltaSNP),
                         REF_FRQ=mean(REF_FRQ),
                         DP.LOW=sum(DP.LOW),
                         DP.HIGH=sum(DP.HIGH)),
                    list(CHROM, POS, REF=REF.x, ALT=ALT.x)]

    m <- m[!is.na(deltaSNP) & deltaSNP!=Inf & deltaSNP!=-Inf]


    df_filt <- filterSNPs(SNPset = as.data.frame(m),
                         refAlleleFreq = 0.05,
                         minTotalDepth = 10,
                         maxTotalDepth = 5000,
                         depthDifference = 1500,
                         minSampleDepth = 10,
                         verbose = TRUE)


     gprime <- runGprimeAnalysis(SNPset=df_filt, windowSize=1000000) #250000
     gprime <- as.data.table(gprime)
     gprime[Gprime>10]

     peaks <- getQTLTable(gprime, alpha=.1)
     peaks <- as.data.table(peaks)
     peaks








    gprime <- as.data.table(gprime)

    summary(gprime[nSNPs>1])

  save(gprime, file="~/gprime_new.Rdata")
  scp aob2x@rivanna.hpc.virginia.edu:~/gprime_new.Rdata ~/.
  R --vanilla
  library(ggplot2); library(data.table)
  load("~/gprime_new.Rdata")
  setnames(gprime, c("CHROM", "POS"), c("chr", "pos"), skip_absent=T)

  ggplot(data=gprime, aes(x=pos, y=Gprime, color=chr)) + geom_point() + facet_grid(~chr)


#### for two separate reps
  setkey(gprime, chr, pos)
  setkey(o, chr, pos)
  o.ag <- o[,list(p.z=mean(p.z)), list(id, chr, pos, set, term, perm)]
  setkey(o.ag, chr, pos)

  m <- merge(gprime, o.ag, allow.cartesian=T)


  overlap <- foreach(perm.i=unique(m$perm), .combine="rbind")%dopar%{
    message(perm.i)
    foreach(term.i=unique(m$term), .combine="rbind")%do%{
      foreach(r=c(1,2), .combine="rbind")%do%{
        #perm.i<-1; term.i="fill"; r<-1
        p1 <- m[perm==perm.i][term==term.i][rep==r]$p.z
        p2 <- m[perm==perm.i][term==term.i][rep==r]$pvalue

        #p1 <- m[perm==perm.i][term==term.i]$p.z
        #p2 <- m[perm==perm.i][term==term.i]$pvalue

        q1 <- rank(p1)/(length(p1)+1)
        q2 <- rank(p2)/(length(p2)+1)


        r#ep1<- fisher.test(table(p1<.01, p2<.01))
        rep1<- fisher.test(table(q1<.05, q2<.05))

        data.table(perm=perm.i, term=term.i, rep=r, or=rep1$estimate, p=rep1$p.value)

        #m[perm==perm.i][term==term.i][rep==r][p.z<.001 & pvalue<.01]

      }
    }
  }

  mean(overlap[rep==2][perm==0][term=="male"]$or > overlap[rep==2][perm!=0][term=="male"]$or)
  mean(overlap[rep==1][perm==0][term=="male"]$or > overlap[rep==1][perm!=0][term=="male"]$or)
  mean(overlap[rep==2][perm==0][term=="fill"]$or > overlap[rep==2][perm!=0][term=="fill"]$or)
  mean(overlap[rep==1][perm==0][term=="fill"]$or > overlap[rep==1][perm!=0][term=="fill"]$or)

  head(overlap[order(p)])
  tail(overlap[order(or)])

  m[,cp:= 1-pchisq(-2*(log(pvalue)+log(p.z)), 4)]
  m[,list(mincp=min(cp)), list(term, perm, rep)][order(mincp)]

  table(m$cp<.000005, m$perm, m$rep)

  save(m, file="~/f1_pool.Rdata")

  scp aob2x@rivanna.hpc.virginia.edu:~/f1_pool.Rdata ~/.

  library(data.table)
  library(ggplot2)

  load("~/f1_pool.Rdata")




#### how many sites
m <- foreach(perm.i=unique(m$perm), .combine="rbind")%do%{
  foreach(term.i=unique(m$term), .combine="rbind")%do%{
    foreach(r=c(1,2), .combine="rbind")%do%{
      #perm.i<-0; term.i="male"; r<-2
      tmp <- m[perm==perm.i][term==term.i][rep==r]

      tmp[,pool.q:=rank(pvalue)/(length(pvalue)+1)]
      tmp[,f1.q:=rank(p.z)/(length(p.z)+1)]
      tmp[,cp:=1-pchisq(-2*(log(pool.q)+log(f1.q)), 4)]
      #m[perm==perm.i][term==term.i][rep==r][p.z<.001 & pvalue<.01]
      tmp[,cp:=pool.q/2 + f1.q/2]

    }
  }
}



  manhattan <- ggplot(data=m[perm==0], aes(x=pos, y=-log10(cp), color=chr)) +
  geom_line() +
  facet_grid(term~chr)


  ggplot(data=m[perm==0], aes(x=pos, y=Gprime, color=chr)) +
  geom_line(data=m[perm==0], aes(x=pos, y=Gprime, color=chr)) +
  geom_point(data=m[perm==0][p.z<.005], aes(x=pos, y=Gprime), color="black") +
  facet_grid(term+rep~chr)








  qtlseq <- runQTLseqAnalysis(SNPset=df_filt, windowSize=250000, bulkSize=40)
  getQTLTable(SNPset=qtlseq, method="QTLseq")

  p <- plotQTLStats(SNPset = gprime, var = "Gprime") +
       geom_hline(aes(yintercept=min(gprime[qvalue<.1]$Gprime)), color="red") +
       geom_hline(aes(yintercept=min(gprime[qvalue<.05]$Gprime)), color="blue") +
       geom_hline(aes(yintercept=min(gprime[qvalue<.01]$Gprime)), color="green")

  #ggsave(p, file="~/D8_QTL.pdf", width=14, height=7)


  peaks <- as.data.table(getQTLTable(SNPset = gprime, method = "Gprime", alpha = 0.05, export = FALSE))
  save(peaks, gprime, file="~/peaks.Rdata")

  #save(gprime, peaks, file="/mnt/sammas_storage/bergland-lab/alan/peaks.Rdata")




### scp aob2x@rivanna.hpc.virginia.edu:~/peaks.Rdata ~/.


library(data.table)
library(ggplot2)

load("~/peaks.Rdata")
gprime

ggplot(data=gprime, aes(x=POS, y=Gprime, color=CHROM))  + geom_line() + facet_grid(~CHROM)

















### LD clumped
### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(ggtern)
  library(viridis)
  library(bigstatsr)
  library(SeqArray)
  library(bigsnpr)
  library(regioneR)


  ### load HWE summary data (HWE_simulations.gather.rivanna.R makes this file)
    load(file="/mnt/sammas_storage/bergland-lab/Daphnia_HWE/hwe.ag.m.ag")

  ### open HWE & genotype dist statistic data
    ### cp /scratch/aob2x/daphnia_hwe_sims/hwe_stat.Rdata /nv/vol186/bergland-lab/Daphnia_HWE/hwe_stat.Rdata
    load(file="/mnt/sammas_storage/bergland-lab/Daphnia_HWE/hwe_stat.Rdata")
    hwe.stat[,py:=paste(pond, year, sep=".")]

  ### load precomputed file (HWE_simulations.prep.rivanna.R makes this file)
    load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
    sc[,year:=tstrsplit(clone, "_")[[2]]]

  ### open genotype file
    genofile <- seqOpen("/mnt/spicy_3/Karen/201620172018FinalMapping/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.seq.gds")

  ### LD blocks
    load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/ldBlocks.Rdata")

  ### test of overlaping ranges
    A <- makeGRangesFromDataFrame(data.frame(chr=ld.blocks[py=="D8.2016.2017.2018.2019"][KB>0]$CHR,
                                        start=ld.blocks[py=="D8.2016.2017.2018.2019"][KB>0]$BP1,
                                        end=ld.blocks[py=="D8.2016.2017.2018.2019"][KB>0]$BP2),
                             start.field="start", end.field="end")

    B <- makeGRangesFromDataFrame(data.frame(chr=peaks$CHROM,
                                       start=peaks$start,
                                       end=peaks$end),
                            start.field="start", end.field="end")

    genome <- getGenomeAndMask(genome=snp.dt[(final.use), list(start=min(pos), end=max(pos)), list(chr)], mask=NA)

    pt <- overlapPermTest(A=A, B=B, genome=genome$genome, ntimes=1000)
    plot(pt)





### nice plots

### fold in haplotype estimates
  ### libraries
    library(data.table)
    library(foreach)
    library(cowplot)
    library(SeqArray)

  ### load output from harp (generated by `plotHarp.R`)
    #load(file="/scratch/aob2x/daphnia_hwe_sims/harp_pools/summarizedOut/harpWins.Rdata")
    load("/mnt/sammas_storage/bergland-lab/alan/harpWins.Rdata")
    setkey(o, chr)
    harpWins <- o
    rm(o)

  ### load sliding window plot (from above)
    load(file="~/slidingWindow_oddsRatio.Rdata")

  ### load precomputed mInform file (from above)
    load(file="~/mInform.Rdata")

  ### open connection to annotated GDS file
    #genofile <- seqOpen("/mnt/pricey_2/DaphniaGenomeAnnotation/snpEff/totalnewmapwMarch2018_Dfiltsnps10bpindels_snps_filter_pass_lowGQmiss.annotated.gds")
    genofile <- seqOpen("/mnt/sammas_storage/bergland-lab/alan/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")
    load("/mnt/sammas_storage/bergland-lab/alan/finalsetsnpset01pulex_20200207.Rdata")

    allSNPs.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                             chr=seqGetData(genofile, "chromosome"),
                             pos=seqGetData(genofile, "position"))


  ### superclone (note, the 'B' superclone is now the 'C' superclone)
    sc <- fread("/mnt/sammas_storage/bergland-lab/alan/CloneInfoFilePulexandObtusa_withmedrd_20200207")

    samps <- seqGetData(genofile, "sample.id")

  ### load haplotype priors
    dat <- fread("/mnt/sammas_storage/bergland-lab/alan/testTrio.consensus.header.phase.csv")

    setnames(dat, c(1,2), c("chr", "pos"))

    dat <- dat[!(A=="1/1" & B=="1/1")][!(A=="0/0" & B=="0/0")]

    dat.p <- dat[,list(ref=REF, alt=ALT,
                        A1=as.numeric(substr(A, 0, 1)),
                        A2=as.numeric(substr(A, 3, 3)),
                        B1=as.numeric(substr(B, 0, 1)),
                        B2=as.numeric(substr(B, 3, 3))),
                    list(chr, pos)]

  ### load pulicaria
    #puli <- fread("/mnt/sammas_storage/bergland-lab/alan/pulicaria.sort.D84a.rmDup.aseReadCounter.allvariant.delim")
    #setnames(puli, c("contig", "position"), c("chr", "pos"))
    #setkey(puli, chr, pos)
    #puli[,puli.geno:=round(refCount/totalCount*10)/10]
    #puli[puli.geno!=1 & puli.geno!=0, puli.geno:=.5]

    #### playing
      #setkey(dat.p, chr, pos)
      #foo <- merge(puli[,c("chr", "pos", "puli.geno")], dat.p, all.x=T, all.y=T)
      #prop.table(table(foo$puli.geno==foo$A.geno))
      #prop.table(table(foo$puli.geno==foo$B.geno))

      #foo.l <- melt(foo, id.vars=c("chr", "pos", "puli.geno", "ref", "alt"))
      #foo.l.ag <- foo.l[puli.geno%in%c(0,1), list(IBS=mean(value==puli.geno, na.rm=T), .N), list(chr, variable)]

      #ggplot(data=foo.l.ag[N>1000], aes(x=variable, y=IBS, color=chr)) + geom_point() + coord_flip()

    load("/mnt/sammas_storage/bergland-lab/alan/puAg.Rdata")
    setkey(pu.ag, chr, pos)



  ### generte windows
    setkey(dat.p, chr, pos)
    chrs <- m.inform[use.chr==T, list(min=min(pos), max=max(pos)), chr]
    window.bp <- 50000
    step.bp <- 5000

    wins <- foreach(i=1:dim(chrs)[1], .combine="rbind")%dopar%{
      #i<-1
      print(i)
      tmp <- data.table(chr=chrs[i]$chr, start=seq(from=chrs[i]$min, to=chrs[i]$max-window.bp, by=step.bp))
      tmp[,stop:=start+window.bp]
      tmp
    }
    wins[,i:=1:dim(wins)[1]]

    setkey(m.inform, chr, pos)


  ### polarized Manhattan plot
    o[,ppa:=(sign(beta.int)*log10(pa))]

  ### plotting functions

    plot.or <- function(i=o[which.min(ppa)]$win.i) {
      #i <- wins[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(start-8660157))]$i
      #i <- o[nSNPs>25][which.max(ppa)]$win.i
      #i <- o[win.i>5000][which.min(ppa)]$win.i

      m.inform[,win:=0]
      m.inform[chr==wins$chr[i] & pos>=wins$start[i] & pos<=wins$stop[i], win:=1]

      ml <- melt(m.inform[,c("chr", "pos", "win", "D8Male.f", "D8PE.f", "BA.geno.diff", "A.geno", "B.geno"), with=F],
                 id.vars=c("chr", "pos", "win", "BA.geno.diff", "A.geno", "B.geno"))
      ml[,winf:=c("genome", "targetWindow")[win+1]]
      ml[,exp_hom:= 1 - 2*value*(1-value)]

      ggplot(m.inform, aes(x=as.factor((BA.geno.diff)), y=log2(or.fold),
                           group=interaction((BA.geno.diff), as.factor(win)), fill=as.factor(win))) +
                 geom_boxplot()
               }

    plot.freq <- function(i=o[which.min(ppa)]$win.i) {
     #i <- wins[chr=="Scaffold_7757_HRSCAF_8726"][which.min(abs(start-8660157))]$i
     #i <- o[nSNPs>25][which.max(ppa)]$win.i
     #i <- o[win.i>5000][which.min(ppa)]$win.i

     m.inform[,win:=0]
     #m.inform[chr==wins$chr[i] & pos>=wins$start[i] & pos<=wins$stop[i], win:=1]
     m.inform[chr==o[win.i==i]$chr & pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop, win:=1]



     ml <- melt(m.inform[,c("chr", "pos", "win", "D8Male.f", "D8PE.f", "BA.geno.diff", "A.geno", "B.geno"), with=F],
                id.vars=c("chr", "pos", "win", "BA.geno.diff", "A.geno", "B.geno"))
     ml[,winf:=c("genome", "targetWindow")[win+1]]
     ml[,exp_hom:= 1 - 2*value*(1-value)]


     ggplot(ml, aes(x=as.factor(variable), y=value,
                           group=interaction(variable, as.factor(winf)), fill=as.factor(winf))) +
                  geom_boxplot() +
                  facet_grid(A.geno~B.geno) +
                  ylab("Alt. allele freq") +
                  theme(axis.text.x = element_text(angle = 45, hjust = 1))
                }

    plot.hap <- function(i=o[which.min(ppa)]$win.i) {

     ### haplotype painting plot
      ### extract out A&B haplotypes
       tmp <- dat.p[chr==o[win.i==i]$chr & pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]

       tmp <- merge(tmp, pu.ag[,c("chr", "pos", "puli.geno"), with=F], all.x=T)
       tmp[,puli.geno:=puli.geno/2]

       tmp.ag <- tmp[,list(i=seq(from=min(pos), to=max(pos), length.out=length(pos)), pos=pos), list(chr)]
       setkey(tmp.ag, chr, pos)
       setkey(tmp, chr, pos)
       tmp <- merge(tmp, tmp.ag)

       tmp[A1==0, A1.allele:=alt]
       tmp[A1==1, A1.allele:=ref]
       tmp[A2==0, A2.allele:=alt]
       tmp[A2==1, A2.allele:=ref]

       tmp[B1==0, B1.allele:=alt]
       tmp[B1==1, B1.allele:=ref]
       tmp[B2==0, B2.allele:=alt]
       tmp[B2==1, B2.allele:=ref]

       tmp[puli.geno==0, puli.allele:=alt]
       tmp[puli.geno==1, puli.allele:=ref]
       tmp[puli.geno==.5, puli.allele:=ref]

       tmp2 <- tmp
       tmp2[tmp$B2==0, A1:=abs(1-A1)]
       tmp2[tmp$B2==0, A2:=abs(1-A2)]
       tmp2[tmp$B2==0, B1:=abs(1-B1)]
       tmp2[tmp$B2==0, B2:=abs(1-B2)]
       tmp2[tmp$B2==0, puli.geno:=abs(1-puli.geno)]



       tmp2[tmp$B2==0 & A1==0, A1.allele:=ref]
       tmp2[tmp$B2==0 & A1==1, A1.allele:=alt]
       tmp2[tmp$B2==0 & A2==0, A2.allele:=ref]
       tmp2[tmp$B2==0 & A2==1, A2.allele:=alt]

       tmp2[tmp$B2==0 & B1==0, B1.allele:=ref]
       tmp2[tmp$B2==0 & B1==1, B1.allele:=alt]
       tmp2[tmp$B2==0 & B2==0, B2.allele:=ref]
       tmp2[tmp$B2==0 & B2==1, B2.allele:=alt]

       tmp2[tmp$B2==0 & puli.geno==0, puli.allele:=ref]
       tmp2[tmp$B2==0 & puli.geno==1, puli.allele:=alt]
       tmp2[tmp$B2==0 & puli.geno==0.5, puli.allele:=ref]



       tmpl <- melt(tmp2, id.vars=c("chr", "pos", "ref", "alt", "i"),
                      measure=list(c("A1", "A2", "B1", "B2", "puli.geno"),
                                   c("A1.allele", "A2.allele", "B1.allele", "B2.allele", "puli.allele")),
                       value.name=c("value", "allele"))
       tmpl[,variable:=c("A1", "A2", "B1", "B2", "puli.geno")[variable]]

     ### fold in annotations

      # THIS line is problematic
       seqSetFilter(genofile, variant.id=allSNPs.dt[chr==o[win.i==i]$chr][pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]$variant.id)


       ann.tmp <- seqGetData(genofile, "annotation/info/ANN")
       ann.tmp.dt <- data.table(times=ann.tmp$length, id=seqGetData(genofile, "variant.id"))
       ann.tmp.l <- ann.tmp.dt[,list(id=rep(id, times)), list(i=id)]
       ann.tmp.l[,ann:=ann.tmp$data]

       classes <- c("stop|mis")
       ann.tmp.l <- ann.tmp.l[grepl(classes, ann)]

       setkey(ann.tmp.l, id)
       setkey(m.inform, id)

       ann.tmp.l[,allele:=tstrsplit(ann, "\\|")[[1]]]
       ann.tmp.l[,mut.class:=tstrsplit(ann, "\\|")[[2]]]
       ann.tmp.l[,gene:=tstrsplit(ann, "\\|")[[4]]]
       ann.tmp.l <- ann.tmp.l[,list(.N), list(id, allele, mut.class, gene)]
       ann.tmp.m <- merge(ann.tmp.l[,c("id","allele", "mut.class")], m.inform[,c("chr", "pos", "id", "A.geno", "B.geno")], by="id")

       setkey(tmpl, chr, pos, allele)
       setkey(ann.tmp.m, chr, pos, allele)

       foo <- merge(tmpl, ann.tmp.m, all.x=T)

       haps <- ggplot() +
                geom_tile(data=foo, aes(x=variable, y=i, fill=as.factor(abs(1-value)))) +
                geom_jitter(data=foo[!is.na(mut.class)], aes(x=variable, y=i, shape=mut.class, size=mut.class), color="black", width=.1) +
                coord_flip()  +
                theme(legend.position = "none")

     ### annotations
      #seqSetFilter(genofile, variant.id=m.inform[chr==o[win.i==i]$chr][pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]$id)
      #ann.tmp <- seqGetData(genofile, "annotation/info/ANN")
      #ann.tmp.dt <- data.table(times=ann.tmp$length, id=seqGetData(genofile, "variant.id"))
      #ann.tmp.l <- ann.tmp.dt[,list(id=rep(id, times)), list(i=id)]
      #ann.tmp.l[,ann:=ann.tmp$data]

      #setkey(ann.tmp.l, id)
      #setkey(m.inform, id)

      #ann.tmp.m <- merge(ann.tmp.l, m.inform[,c("chr", "pos", "id", "A.geno", "B.geno")])

      #ann.tmp.m[grepl("stop", ann)]$pos[1]
      #ann.tmp.m[grepl("missense", ann)]

     ### haplotype frequency plot
      pad <- 0
      o[,mid:=start/2 + stop/2]
      hapFreq.tmp <- harpWins

      hapFreq.tmp[,targetWin:=0]
      hapFreq.tmp <- hapFreq.tmp[chr==o[win.i==i]$chr & start<=(o[win.i==i]$mid-pad) & stop>=(o[win.i==i]$mid+pad), targetWin:=1]

      hapFreq.tmp.l <- melt(hapFreq.tmp, id.vars=c("chr", "start", "stop", "pool", "targetWin"))
      hapFreq.tmp.l[grepl("Male", pool),pool.f:="D8Male.f"]
      hapFreq.tmp.l[grepl("PE", pool),pool.f:="D8PE.f"]

      hapFreq.tmp.lag.ag <- hapFreq.tmp.l[,list(freq=mean(value, na.rm=T), sd=sd(value, na.rm=T)), list(allele=variable, pool=factor(pool, levels=rev(c("D8PE1", "D8PE2", "D8Male1", "D8Male2"))), targetWin)]

      haplofreq.plot <- ggplot() +
      geom_bar(data=hapFreq.tmp.lag.ag[targetWin==1], aes(x=allele, y=freq, group=pool, fill=pool), stat="identity", position="dodge") +
      geom_point(data=hapFreq.tmp.lag.ag[targetWin==0], aes(x=allele, y=freq, group=pool), position=position_dodge(1)) +
      geom_errorbar(data=hapFreq.tmp.lag.ag[targetWin==0], aes(x=allele, ymin=freq-sd, ymax=freq+sd, group=pool), position=position_dodge(1)) +
      coord_flip() + guides(fill = guide_legend(reverse = TRUE))

      ### map
        map <- ggplot() +
        geom_point(data=tmp, aes(y=1, x=seq(from=min(pos), to=max(pos), length.out=length(pos)))) +
        geom_point(data=tmp, aes(y=0, x=pos)) +
        geom_segment(data=tmp, aes(y=0, yend=1, x=pos, xend=seq(from=min(pos), to=max(pos), length.out=length(pos))), size=.01) +
        scale_size_manual()


      ### density
        dens <- ggplot() +
                geom_histogram(data=o, aes(nSNPs)) +
                geom_vline(data=o[win.i==i], aes(xintercept=nSNPs), color="red")


     plot_grid(haps, haplofreq.plot, map, dens, ncol=2)

    }

    plot.hap.freq <- function(i=o[which.min(ppa)]$win.i, bi, aj) {
      message(paste(i, bi, aj, sep=" / "))
      #i <- o[which.min(ppa)]$win.i
      #
      #x11()
      #bi<-1; aj<-0


      #rm(pooled.tmp, haps.tmp, pool.hap)

      pooled.exp <- m.inform[,list(freq.mu=c(mean(D8Male.f, na.rm=T), mean(D8PE.f, na.rm=T)), variable=c("D8Male.f", "D8PE.f")), list(A.geno, B.geno)]

      pooled.tmp <- m.inform[chr==o[win.i==i]$chr & pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]
      haps.tmp <- dat.p[chr==o[win.i==i]$chr & pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]


      setkey(pooled.tmp, chr, pos)
      setkey(haps.tmp, chr, pos)
      pool.hap <- merge(pooled.tmp[,c("chr", "pos", "A.geno", "B.geno", "D8Male.f", "D8PE.f"), with=F], haps.tmp, all=T)

      pool.hap[,A.haploSum:=abs(2-(A1+A2))]
      pool.hap[,B.haploSum:=abs(2-(B1+B2))]

      print(table(A=pool.hap$A.geno, B=pool.hap$B.geno))

      ### check
        #table(pool.hap$A.geno==pool.hap$A.haploSum); table(pool.hap$B.geno==pool.hap$B.haploSum)

      Bi.Aj <- na.omit(pool.hap[B.geno==bi][A.geno==aj])
      message(dim(Bi.Aj)[1])

      if(dim(Bi.Aj)[1]!=0) {
        Bi.Aj.freq.l <- melt(Bi.Aj, id.vars=c("chr", "pos"), measure.vars=c("D8Male.f", "D8PE.f"))
        Bi.Aj.haps.l <- melt(Bi.Aj, id.vars=c("chr", "pos"), measure.vars=c("A1", "A2", "B1", "B2"))

        freq.plot <- ggplot() +
                     geom_boxplot(data=Bi.Aj.freq.l, aes(x=as.factor(variable), y=value,
                                           group=variable)) +
                     ylab("Alt. allele freq") +
                     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                     ggtitle(paste("A:", aj, " B:", bi, "   nSNPs: ", dim(Bi.Aj)[1], sep="")) +
                     ylim(0,1) +
                     geom_point(data=pooled.exp[A.geno==aj & B.geno==bi], aes(x=as.factor(variable), y=freq.mu), color="red")

        hap.plot <- ggplot(data=Bi.Aj.haps.l, aes(x=variable, y=i, fill=as.factor(abs(1-value)))) + geom_tile() + coord_flip()  + theme(legend.position = "top")
        return(plot_grid(freq.plot, hap.plot))
      } else {
        freq.plot <- ggplot(data=data.table(x=1, y=1), aes(x=x, y=x))
        hap.plot <- freq.plot
        return(plot_grid(freq.plot, hap.plot))
      }

    }

    plot.genos <- function(i=o[win.i<5000][which.min(ppa)]$win.i, ii) {

      #i=o[win.i<5000][which.min(ppa)]$win.i



      region <- m.inform[chr==o[win.i==i]$chr][pos>=o[win.i==i]$start & pos<=o[win.i==i]$stop]
      setkey(region, chr, pos)
      setkey(allSNPs.dt, chr, pos)

      sets <- list(A_v_C=sc[SC%in%c("A", "C")]$clone,
                   D8=sc[population=="D8"]$clone,
                   DBunk=sc[population=="DBunk"]$clone)



      seqResetFilter(genofile)
      seqSetFilter(genofile,
                   sample.id=sets[[ii]],
                   variant.id=allSNPs.dt[J(region)]$variant.id)

      geno.mat <- as.data.table(seqGetData(genofile, "$dosage"))
      setnames(geno.mat,
                names(geno.mat),
               paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), c(1:length(seqGetData(genofile, "position"))), sep=";"))
      geno.mat[,clone:=seqGetData(genofile, "sample.id")]
      gml <- melt(geno.mat, id.vars="clone")
      gml[,chr:=tstrsplit(variable, ";")[[1]]]
      gml[,pos:=tstrsplit(variable, ";")[[2]]]
      gml[,i:=tstrsplit(variable, ";")[[3]]]
      setnames(gml, "value", "geno")

      setkey(gml, "clone")
      setkey(sc, "clone")
      gml <- merge(gml, sc)
      gml.ag <- gml[,list(r=superClone.sizeRank[1], clone=unique(clone)), list(SC)]

      gml[,clonef:=factor(clone, levels=gml.ag[order(r)]$clone)]

      #if(ii==1) gml[,clonef:=factor(clone, levels=sc[SC%in%c("A", "C")][order(SC)]$clone)]
      #if(ii==2) gml[,clonef:=factor(clone, levels=sc[population=="D8"][order(SC)]$clone)]


      gml[,clonefn:=as.numeric(clonef)]

      ggplot(data=gml, aes(y=i, x=clonefn, fill=as.factor(geno))) + geom_tile() + coord_flip() + facet_wrap(~SC, scales="free")
    }


  ### QTL peaks
    load(file="/mnt/sammas_storage/bergland-lab/alan/peaks.Rdata")

    i <- 6
    plot.hap(o[chr==peaks[i]$CHROM][which.min(sqrt((peaks[i]$start/2 + peaks[i]$end/2 -start)^2 + (peaks[i]$start/2 + peaks[i]$end/2 -stop)^2))]$win.i)
















### fold in HWE stuff
  ### libraries
    library(SeqArray)
    library(data.table)

  ### funciton to radomly subset one clone per superclone
    subsampClone <- function(sc.dt, n=1, use.pond="DBunk") {
      sc.samp <- sc.dt[,list(clone=sample(clone, size=n)), list(sc.uniq)]
      sc.samp[,pond:=tstrsplit(clone, "_")[[3]]]
      sc.samp[,year:=tstrsplit(clone, "_")[[2]]]
      return(sc.samp[pond%in%use.pond])
    }


  ### load QTL peaks
    load(file="/mnt/sammas_storage/bergland-lab/alan/peaks.Rdata")

  ### load older phased VCF file that
    phased <- seqOpen("/mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.biallelic.phased.gds")

  ### load superclone data
    load("/mnt/spicy_3/AlanDaphnia/outputData/superclone.Rdata")
    sc <- fread("/mnt/spicy_3/AlanDaphnia/vcf/Superclones20161718updated_20190802.csv")
    sc[,pond := tstrsplit(clone, "_")[[3]]]
    sc[,sc.uniq := SC]
    sc[SC=="OO", sc.uniq:=paste(SC, SCnum, sep=".")]

  ### snp.dt
    snp.dt <- data.table(chr=seqGetData(phased, "chromosome"),
                         pos=seqGetData(phased, "position"),
                         id=seqGetData(phased, "variant.id"))


   ### test
    i <- which.min(peaks$meanQval)
    sub.sc <- subsampClone(sc.dt=sc, n=1, use.pond="D8")

    seqResetFilter(phased)
    seqSetFilter(phased,
                 variant.id=snp.dt[chr==peaks[i]$CHROM][pos>=peaks[i]$start][pos<=peaks[i]$end]$id)

    seqGetData(phased, "phase")












### VarScan data
### libraries
  library(data.table)
  library(QTLseqr)
  library(ggplot2)

### pooled data
  load("/mnt/sammas_storage/bergland-lab/alan/totalADRDlongall.Rdata")

  geno[,effRD:=floor(RRD)]
  geno[,effPA:=round(propalt*effRD)/effRD]


  #### D8
    geno.w <- dcast(geno[pond=="DBunk"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

    geno.w.ag <- geno.w[,.N,chr]
    geno.w <- geno.w[chr%in%geno.w.ag[N>1000]$chr]

    geno.w.comb <- geno.w[,list(AD_REF.LOW = (1-effPA_D8Male1)*effRD_D8Male1 + (1-effPA_D8Male2)*effRD_D8Male2,
                                AD_ALT.LOW = (effPA_D8Male1)*effRD_D8Male1 + (effPA_D8Male2)*effRD_D8Male2,
                                AD_REF.HIGH = (1-effPA_D8PE1)*effRD_D8PE1 + (1-effPA_D8PE2)*effRD_D8PE2,
                                AD_ALT.HIGH = (effPA_D8PE1)*effRD_D8PE1 + (effPA_D8PE2)*effRD_D8PE2),
                            list(CHROM=chr, POS=pos)]

    geno.w.comb[,deltaSNP:=(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - (AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
    geno.w.comb[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
    geno.w.comb[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
    geno.w.comb[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]

    geno.w.comb.filt <- filterSNPs(SNPset = as.data.frame(geno.w.comb),
                         refAlleleFreq = 0.05,
                         minTotalDepth = 50,
                         maxTotalDepth = 350,
                         depthDifference = 200,
                         minSampleDepth = 25,
                         verbose = TRUE)


    gprime <- runGprimeAnalysis(SNPset=geno.w.comb.filt, windowSize=250000, outlierFilter="Hampel", maxk=100)
    gprime <- as.data.table(gprime)
    p <- plotQTLStats(SNPset = gprime, var = "Gprime") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.1]$Gprime)), color="red") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.05]$Gprime)), color="blue") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.01]$Gprime)), color="green") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.001]$Gprime)), color="green")

    p
    peaks <- as.data.table(getQTLTable(SNPset = gprime, method = "Gprime", alpha = 0.0001, export = FALSE))
    peaks[order(meanGprime)]


  ### DBunk
    geno.w <- dcast(geno[pond=="DBunk"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

    geno.w.ag <- geno.w[,.N,chr]
    geno.w <- geno.w[chr%in%geno.w.ag[N>1000]$chr]

    geno.w.comb <- geno.w[,list(AD_REF.LOW = (1-effPA_DBunkMale)*effRD_DBunkMale,
                                AD_ALT.LOW = (effPA_DBunkMale)*effRD_DBunkMale,
                                AD_REF.HIGH = (1-effPA_DBunkPE1)*effRD_DBunkPE1 + (1-effPA_DBunkPE2)*effRD_DBunkPE2,
                                AD_ALT.HIGH = (effPA_DBunkPE1)*effRD_DBunkPE1 + (effPA_DBunkPE2)*effRD_DBunkPE2),
                            list(CHROM=chr, POS=pos)]

    geno.w.comb[,deltaSNP:=(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - (AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
    geno.w.comb[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
    geno.w.comb[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
    geno.w.comb[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]



    geno.w.comb.filt <- filterSNPs(SNPset = as.data.frame(geno.w.comb),
                         refAlleleFreq = 0.05,
                         minTotalDepth = 50,
                         maxTotalDepth = 350,
                         depthDifference = 200,
                         minSampleDepth = 25,
                         verbose = TRUE)


    gprime <- runGprimeAnalysis(SNPset=geno.w.comb.filt, windowSize=250000, outlierFilter="Hampel", maxk=100)
    gprime <- as.data.table(gprime)
    p <- plotQTLStats(SNPset = gprime, var = "Gprime") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.1]$Gprime)), color="red") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.05]$Gprime)), color="blue") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.01]$Gprime)), color="green") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.001]$Gprime)), color="green")

    p
    peaks <- as.data.table(getQTLTable(SNPset = gprime, method = "Gprime", alpha = 0.0001, export = FALSE))
    peaks[order(meanGprime)]

  ### DCat
    geno.w <- dcast(geno[pond=="DCat"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

    geno.w.ag <- geno.w[,.N,chr]
    geno.w <- geno.w[chr%in%geno.w.ag[N>1000]$chr]

    geno.w.comb <- geno.w[,list(AD_REF.LOW = (1-effPA_DCatMale)*effRD_DCatMale,
                                AD_ALT.LOW = (effPA_DCatMale)*effRD_DCatMale,
                                AD_REF.HIGH = (1-effPA_DCatPE1)*effRD_DCatPE1 + (1-effPA_DCatPE2)*effRD_DCatPE2,
                                AD_ALT.HIGH = (effPA_DCatPE1)*effRD_DCatPE1 + (effPA_DCatPE2)*effRD_DCatPE2),
                            list(CHROM=chr, POS=pos)]

    geno.w.comb[,deltaSNP:=(AD_ALT.LOW/(AD_ALT.LOW+AD_REF.LOW)) - (AD_ALT.HIGH/(AD_ALT.HIGH+AD_REF.HIGH))]
    geno.w.comb[,REF_FRQ:=(AD_REF.LOW + AD_REF.HIGH) / (AD_REF.LOW + AD_REF.HIGH + AD_ALT.LOW + AD_ALT.HIGH )]
    geno.w.comb[,DP.LOW:=AD_REF.LOW + AD_ALT.LOW]
    geno.w.comb[,DP.HIGH:=AD_REF.HIGH + AD_ALT.HIGH]



    geno.w.comb.filt <- filterSNPs(SNPset = as.data.frame(geno.w.comb),
                         refAlleleFreq = 0.05,
                         minTotalDepth = 50,
                         maxTotalDepth = 350,
                         depthDifference = 200,
                         minSampleDepth = 25,
                         verbose = TRUE)



    gprime <- runGprimeAnalysis(SNPset=geno.w.comb.filt, windowSize=250000, outlierFilter="Hampel", maxk=100)
    gprime <- as.data.table(gprime)
    p <- plotQTLStats(SNPset = gprime, var = "Gprime") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.1]$Gprime)), color="red") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.05]$Gprime)), color="blue") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.01]$Gprime)), color="green") +
         geom_hline(aes(yintercept=min(gprime[qvalue<.001]$Gprime)), color="green")

    p
    peaks <- as.data.table(getQTLTable(SNPset = gprime, method = "Gprime", alpha = 0.0001, export = FALSE))
    peaks[order(meanGprime)]






























ggplot(data=gprime, aes(x=POS, y=Gprime)) + geom_line() + facet_wrap(~CHROM, nrow=1)

getQTLTable(SNPset = gprime, method = "Gprime",alpha = 0.005, export = FALSE)


gprime <- as.data.table(gprime)
gprime[,ed4:=sqrt(deltaSNP^2 + (-1*deltaSNP)^2)^4]
ggplot(data=gprime, aes(x=POS, y=ed4)) + geom_line() + facet_wrap(~CHROM, nrow=1)

o <- runQTLseqAnalysis(df_filt, windowSize = 500000,
popStruc = "F2",
bulkSize = c(70, 100), replications = 10000,
intervals = c(95, 99) )

plotQTLStats(SNPset = o, var = "deltaSNP", plotIntervals = TRUE))


 p3 <- plotQTLStats(SNPset = gprime, var = "Gprime", plotThreshold = TRUE, q = 0.1)
 p3


 gprime <- as.data.table(gprime)
 gprime[qvalue<.01]
