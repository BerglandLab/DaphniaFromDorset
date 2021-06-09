### scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/justPheno.rQTL.csv ~/.


### libraries
  library(data.table)
  library(qtl)
  library(ggplot2)
  library(cowplot); theme_set(theme_cowplot())
  library(patchwork)

### load peaks
  load("~/peaks.Rdata")
  setnames(peaks, "CHROM", "chr")
  setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))

### load cross
  AxCF1 <- read.cross("csv","",
                      file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/justPheno.rQTL.csv",
                      crosstype="4way", genotypes=NULL)

  AxCF1 <- read.cross("csv","",
                      file="~/justPheno.rQTL.csv",
                      crosstype="4way", genotypes=NULL)

### old  vs new phenotypes
  load("~/mphenoupdate_20200824.Rdata")
  pheno <- as.data.table(AxCF1$pheno)

  ph2 <- merge(mphenoupdate, pheno, by="SCB")

  plot(embresidB~fill.ranef, ph2)
  plot(malesresid~propmalenoneo.ranef, ph2)


### make/fit qtl using peaks from PoolSeq
  AxCF1 <- sim.geno(AxCF1, n.draws=5)

  mq <- makeqtl(AxCF1, chr=peaks$chr, pos=peaks$posMaxGprime, qtl.name=paste("QTL", c(1:dim(peaks)[1]), sep=""))
  mq <- makeqtl(AxCF1, chr=max(mr.fill)$chr, pos=max(mr.fill)$pos)

  lod <- fitqtl(AxCF1, pheno.col=4, mq,
                formula=paste("y~", paste(paste("Q", c(1:12), sep=""), collapse="+ "), sep=" "),
                method="imp")
  summary(lod)


### marker regressoion
  mr.propmale <- scanone(AxCF1, pheno.col=3, method="imp")
  mr.fill <- scanone(AxCF1, pheno.col=4, method="imp")
  mr.epp <- scanone(AxCF1, pheno.col=5, method="imp")
  mr.propmale.blup <- scanone(AxCF1, pheno.col=6, method="mr")
  mr.fill.blup <- scanone(AxCF1, pheno.col=7, method="mr")
  mr.epp.blup <- scanone(AxCF1, pheno.col=8, method="mr")



### mega plot
mr.epp.dt <- as.data.table(mr.epp)
mr.fill.dt <- as.data.table(mr.fill)
mr.propmale.dt <- as.data.table(mr.propmale)

eppresid.plot <- ggplot() +
#geom_hline(yintercept=summary(perm)[2]) +
geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
geom_rect(data=peaks, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="grey", alpha=.5) +
geom_line(data=mr.epp.dt, aes(x=pos, y=lod, color=chr), size=1) +
facet_grid(.~chr, scales="free_x") +
theme(strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
ggtitle("Num. Epp")


#setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))
Gprime.plot <- ggplot(data=gprime, aes(x=pos, y=Gprime, color=chr)) +
#geom_rect(data=peaks, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="grey", alpha=.5) +
geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
geom_line(size=.75) +
facet_grid(.~chr, scales="free_x") +
theme(strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
ggtitle("pooledWild")

embresid.plot <- ggplot() +
#geom_hline(yintercept=summary(perm)[2]) +
geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
geom_rect(data=peaks, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="grey", alpha=.5) +
geom_line(data=mr.fill.dt, aes(x=pos, y=lod, color=chr), size=1) +
facet_grid(.~chr, scales="free_x") +
theme(strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
ggtitle("Fill Rate")


propmale.plot <- ggplot() +
#geom_hline(yintercept=summary(perm)[2]) +
geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
geom_rect(data=peaks, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill="grey", alpha=.5) +
geom_line(data=mr.propmale.dt, aes(x=pos, y=lod, color=chr), size=1) +
facet_grid(.~chr, scales="free_x") +
theme(strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
ggtitle("Prop Male")



### final plot
Gprime.plot / propmale.plot / eppresid.plot / embresid.plot










### Genotype frequencies
  gt1 <- as.data.table(geno.table(AxCF1, scanone.output=T))
  gt2 <- as.data.table(geno.table(AxCF1, scanone.output=F))


  gt <- cbind(gt1[,c("chr", "pos"), with=F], gt2[,c("AC", "BC", "AD", "BD", "P.value"), with=F])
  gt[,n:=AC+BC+AD+BD]
  gt[,AC.exp:=.25*n]
  gt[,BC.exp:=.25*n]
  gt[,AD.exp:=.25*n]
  gt[,BD.exp:=.25*n]

  gt[,chisq:=((AC-AC.exp)^2)/AC.exp +
             ((BC-BC.exp)^2)/BC.exp +
             ((AD-AD.exp)^2)/AD.exp +
             ((BD-BD.exp)^2)/BD.exp]

  gt[,p:=1-pchisq(chisq, 3)]

  gt[chr==max(mr.fill.blup)$chr & pos==max(mr.fill.blup)$pos]

  ggplot(data=gt, aes(x=pos, y=-log10(p))) + geom_point() + facet_wrap(~chr)






  ### QC
    #plot(AxCF1)

  ### reduce to grid
   #newmap <- est.map(AxCF1)
   #plotMap(AxCF1, newmap)
   #AxCF1 <- replace.map(AxCF1, newmap)

   #AxCF1 <- est.rf(AxCF1)
   #plotRF(AxCF1)

   #AxCF1 <- calc.genoprob(AxCF1, stepwidth="fixed")
   #AxCF1.sub <- reduce2grid(AxCF1)
   #AxCF1 <- sim.geno(AxCF1)

  ### marker regression
    ### average
      mr.propmale <- scanone(AxCF1, pheno.col=3, method="mr")
      mr.propmale.perm <- scanone(AxCF1, pheno.col=3, method="mr", n.perm=100)
      plot(mr.propmale)
      abline(h=summary(mr.propmale.perm)[1,1])

      mr.propmale <- scanone(AxCF1, pheno.col=3, method="mr")
      mr.fill <- scanone(AxCF1, pheno.col=4, method="mr")
      mr.fill.perm <- scanone(AxCF1, pheno.col=4, method="mr", n.perm=100)
      #mr.fill.s2<- scantwo(AxCF1, pheno.col=4, method="mr")

      plot(mr.fill)
      abline(h=summary(mr.fill.perm)[1,1])



      plotPXG(AxCF1, pheno.col=4, marker=find.marker(AxCF1,
              peaks[which.max(maxGprime)]$CHROM,
              peaks[which.max(maxGprime)]$posMaxGprime))


                  mr.propmale <- scanone(AxCF1, pheno.col=3, method="mr")
                  mr.fill <- scanone(AxCF1, pheno.col=4, method="mr")
      mr.epp <- scanone(AxCF1, pheno.col=5, method="mr")
      plot(mr.epp)







          plot(mr.propmale)

          effectplot(AxCF1,
                      pheno.col=6,
                      mname1=find.marker(AxCF1,
                                  max(mr.propmale.blup)$chr,
                                  max(mr.propmale.blup)$pos))

          plotPXG(AxCF1,
                      pheno.col=6,
                      marker=find.marker(AxCF1,
                                  max(mr.propmale.blup)$chr,
                                  max(mr.propmale.blup)$pos))



          mr.fill.blup.perm <- scanone(AxCF1, pheno.col=7, method="mr", n.perm=100)

          plot(mr.fill.blup)
          abline(h=summary(mr.fill.blup.perm)[1,1])
          plotPXG(AxCF1,
                      pheno.col=7,
                      marker=find.marker(AxCF1,
                                  max(mr.fill.blup)$chr,
                                  max(mr.fill.blup)$pos))



          mr.epp.blup <- scanone(AxCF1, pheno.col=8, method="mr")
          mr.epp.blup.perm <- scanone(AxCF1, pheno.col=8, method="mr", n.perm=100)

          plot(mr.epp.blup)
          abline(h=summary(mr.epp.blup.perm)[1,1])
