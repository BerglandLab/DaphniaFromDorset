# module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


args = commandArgs(trailingOnly=TRUE)
permUse <- as.numeric(args[1]) - 1
crossType <- args[2]

#permUse <-0; crossType<-"CxC"
### libraries

  library(data.table)
  library("Ckmeans.1d.dp")
  library(foreach)
  library(regioneR)
  library(doMC)
  registerDoMC(10)

################
### computer ###
################

  #  ### setwd
  #    setwd("/Users/alanbergland/Documents/GitHub")
  #
  #  ### load F1 cross data (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`)
  #    load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")
  #
  #  ### load in PoolSeq data (made by `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLSeqR.R`)
  #    load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/gprime_peaks.replicates.250K.05.Rdata")
  #
  #  ### load in F1 mapping data (made by `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/lme4qtl.collect.R`)
  #    load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")

###############
### rivanna ###
###############

  ### set wd
    setwd("/scratch/aob2x/daphnia_hwe_sims/")

  ### load F1 cross data (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`)
    load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")

  ### load in PoolSeq data (made by `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLSeqR.R`)
    load("/project/berglandlab/alan/gprime_peaks.replicates.250K.05.Rdata")

  ### load in F1 mapping data
    load(paste("/project/berglandlab/alan/lme4qtl_output.", crossType, ".long.Rdata", sep=""))


### some data prep for the F1 stuff
  ### we are takign the average just because of a consequence of how the data are output. e.g., `o[id==111][term=="male"][perm==0]`
  o.ag <- o[, list(p.z=mean(p.z, na.rm=T), chisq=mean(chisq), p.aov=mean(p.aov, na.rm=T)), list(term, perm, id, chr, pos)]
  #o.ag.ag <- o.ag[,list(pr=1-mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id, chr, pos)]
  #o.ag.ag[pr==0, pr:=1/201]

  ### define permutation thresholds
  o.ag.perm <- o.ag[perm!=0, list(p.aov.thr=c(quantile(p.aov, .01), quantile(p.aov, .05)),
                                  q=c(.01, .05)), list(term, chr)]

  setkey(o.ag.perm, term, chr)
  #setkey(o.ag, term, chr)

### make genome object
  mask.bed <- fread("/project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.bed")
  setkey(mask.bed, "V1")

  mask.bed <- mask.bed[J(unique(o.ag.perm$chr))]
  genome <- getGenomeAndMask(genome=gprime[, list(start=min(POS), end=max(POS)), list(chr=CHROM)],
                            mask=mask.bed)


### iterate through permutations
  setkey(o.ag, perm)

  overlap.out <- foreach(perm.i=permUse, .combine="rbind", .errorhandling="remove")%dopar%{
    message(perm.i)
    ### get subset
      o.ag.temp <- o.ag[J(perm.i)]
      setkey(o.ag.temp, term, chr)

      o.ag.temp <- merge(o.ag.temp, o.ag.perm[q==0.01], allow.cartesian=T)

    ### identify QTL boundaries from F1 mapping data
      o.ag.id <- foreach(chr.i=unique(o.ag.temp$chr), .combine="rbind")%do%{
        #chr.i<-"Scaffold_2158_HRSCAF_2565"
        tmp <- o.ag.temp[chr==chr.i][,list(pos=unique(pos)), list(chr)]
        setkey(tmp, chr, pos)
        tmp[,id:=c(1:dim(tmp)[1])]
      }
      setkey(o.ag.id, chr, pos)
      setkey(o.ag.temp, chr, pos)
      o.ag.temp <- merge(o.ag.temp, o.ag.id)

      clusters <- o.ag.temp[,
                            list(cluster=Cksegs.1d.dp(id.y[p.aov<=p.aov.thr])$cluster,
                                id.y=id.y[p.aov<=p.aov.thr],
                                pos=pos[p.aov<=p.aov.thr]),
                             list(term, chr)]

      cluster.ag <- clusters[,list(start=min(pos), end=max(pos)), list(term, chr, cluster)]


    ### test of overlaping ranges
      overlap.o <- foreach(pheno=c("fill", "male"), .combine="rbind", .errorhandling="remove")%do%{
        #pheno <- "male"
        A <- makeGRangesFromDataFrame(data.frame(chr=cluster.ag[term==pheno]$chr,
                                            start=cluster.ag[term==pheno]$start,
                                            end=cluster.ag[term==pheno]$end),
                                 start.field="start", end.field="end")


        B <- makeGRangesFromDataFrame(data.frame(chr=peaks$CHROM,
                                           start=peaks$start,
                                           end=peaks$end),
                                start.field="start", end.field="end")

        #suppressWarnings()
        pt <- overlapPermTest(A=A, B=B, genome=genome$genome, ntimes=1000)
        data.table(pheno=pheno, p=pt[[1]]$pval, z=pt[[1]]$zscore, perm=perm.i)
      }

    ### return
        return(overlap.o)
    }


  write.csv(overlap.out, file=paste("/scratch/aob2x/daphnia_hwe_sims/overlap_test/overlap_out_", permUse, ".", crossType, ".csv", sep=""))
