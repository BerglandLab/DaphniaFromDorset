### ijob -c1 -p standard -A berglandlab
### module load intel/18.0 intelmpi/18.0 R/3.6.0; R

### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(foreach)
  library(doMC)
  registerDoMC(4)

### load QTLSeqR output
  #load("/mnt/sammas_storage/bergland-lab/alan/harp_summarized_play.Rdata")
  load(file="/mnt/sammas_storage/bergland-lab/alan/peaks.Rdata")
  setnames(peaks, "CHROM", "chr")

### load pre-computed file (probably from here: DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/harp_windows.R)
  #load("/scratch/aob2x/daphnia_hwe_sims/harp_pools/summarizedOut/harpWins_250K.Rdata")
  load("/mnt/sammas_storage/bergland-lab/alan/harpWins_250K.Rdata")

### a bit of summary of the harp output
  o[,sex:=gsub("1|2", "", pool)]

  o.ag <- o[,list(freq=c(mean(pA1), mean(pA2), mean(pB1), mean(pB2)),
                  allele=c("A1", "A2", "B1", "B2")),
             list(chr, start, stop, sex)]

  o.ag.ag <- o.ag[,list(delta=(freq[sex=="D8Male"]) - (freq[sex=="D8PE"])),
                   list(chr, start, stop, allele)]

  m <- o.ag.ag[,list(dist=sum(delta^2), eff=sum(delta[grepl("A", allele)]) - sum(delta[grepl("B", allele)])), list(chr, start, stop)]

  m[,i:=c(1:dim(m)[1])]

  setkey(m, chr, start, stop)
  setkey(o.ag, chr, start, stop)
  setkey(o.ag.ag, chr, start, stop)

  summary(m)


### merge QTL & peaks files into one plot
  ### First, identify position i position of peaks
   findPeak <- function(i) {
      #i<-1
      setkey(m, chr)
      tmp <- as.data.frame(m)
      tmp <- as.data.table(tmp[m$chr==peaks$chr[i],])
      tmp$posDist <- sqrt((tmp$start - peaks[i]$posMaxGprime)^2 + (tmp$stop - peaks[i]$posMaxGprime)^2)
      which(peaks[i]$posMaxGprime==tmp[which.min(posDist)]$start:tmp[which.min(posDist)]$stop) / (tmp[which.min(posDist)]$stop - tmp[which.min(posDist)]$start) + 0.5 + tmp[which.min(posDist)]$i

    }

    peaks[,i:=unlist(sapply(c(1:dim(peaks)[1]), findPeak))]

    write.csv(peaks, "~/peaks.csv", quote=F, row.names=F)

  ### determine i for each position in Gprime plot
    setkey(gprime, CHROM, POS)
    gprime.i <- foreach(i=m$i)%dopar%{
      print(paste(i, max(m$i)))
      tmp <- gprime[J(data.table(CHROM=m$chr[i], POS=m$start[i]:m$stop[i])), nomatch=0]
      tmp[,i:=(seq(from=(i-1), to=i, length.out=dim(tmp)[1]) + 0.5) ]
      return(tmp)
    }
    gprime.i <- rbindlist(gprime.i)
    setnames(gprime.i, "CHROM", "chr")


  #### individuals QTL effects
    qtl.eff <- foreach(i=1:dim(peaks)[1])%do%{
      buffer <- 50000
      #i<-6

      tmp <- o.ag.ag[chr==peaks$chr[i]]
      tmp[,x:=c(1:dim(tmp)[1])]

      start.x = tmp[which.min(abs(start - (peaks$start[i]-buffer)))]$x
      end.x = tmp[which.min(abs(stop - (peaks$end[i]+buffer)))]$x


      sort(c(start.x:end.x))


      tmp.ag <- tmp[sort(c(start.x:end.x)),list(delta.mu=median(delta), delta.sd=sd(delta), qtl=i, qtl.chr=peaks[i]$chr), list(allele)]
      tmp.ag
    }
    qtl.eff <- rbindlist(qtl.eff)
    qtl.eff[,parent:=substr(allele, 0, 1)]
    qtl.eff.ag <- qtl.eff[,list(delta.mu=sum(delta.mu), delta.sd=sum(delta.sd)), list(qtl, allele=parent)]

  ### combined

    qtl.plot <- ggplot(data=gprime.i, aes(x=i, y=Gprime, color=chr)) +
    geom_vline(data=peaks, aes(xintercept=i), color="black", linetype="dashed") +
    geom_line(size=.75) +
    theme(legend.position = "none")

    effect.plot <- ggplot(data=m, aes(x=i, y=-100*dist, color=chr)) +
    geom_vline(data=peaks, aes(xintercept=i), color="black", linetype="dashed") +
    geom_line(size=.75) +
    theme(legend.position = "none") +
    scale_x_continuous(position = "top")

    allele.plot <- ggplot() +
    geom_bar(data=qtl.eff, aes(x=allele, y=delta.mu, group=allele, fill=allele), stat="identity", position="dodge") +
    facet_grid(~qtl) +
    theme(panel.spacing = unit(0.4, "lines"),
          panel.background=element_rect(fill="lightgrey"),
          panel.border=element_rect(colour="black",size=1))


    locus.sum.plot <- ggplot() +
    geom_bar(data=qtl.eff.ag, aes(x=allele, y=delta.mu, group=allele, fill=allele), stat="identity", position="dodge") +
    facet_grid(~qtl) +
    theme(panel.spacing = unit(0.4, "lines"),
          panel.background=element_rect(fill="lightgrey"),
          panel.border=element_rect(colour="black",size=1))


    big.plot <- plot_grid(qtl.plot, effect.plot, allele.plot, locus.plot, labels=c("A", "B", "C", "D"), nrow=4)
    ggsave(big.plot, file="~/big_plot.png", width=18, height=12)









           plotPeak <- function(i, buffer=5000) {
            #i<-6; buffer<-5000
             tmp1 <- o.ag[chr==peaks[i]$CHROM][start>=(peaks[i]$start-buffer)][stop<=(peaks[i]$end+buffer)]
             tmp2 <- o.ag.ag[chr==peaks[i]$CHROM][start>=(peaks[i]$start-buffer)][stop<=(peaks[i]$end+buffer)]
             tmp3 <- m[chr==peaks[i]$CHROM][start>=(peaks[i]$start-buffer)][stop<=(peaks[i]$end+buffer)]

             deltaPlot <- ggplot(data=tmp2, aes(x=(start/2+stop/2), y=delta, group=allele, color=allele, linetype=as.factor(grepl("A", allele)))) + geom_line()
             freqPlot <- ggplot(data=tmp1, aes(x=(start/2+stop/2), y=freq, group=interaction(allele, sex), color=sex, linetype=allele)) + geom_line()
             distPlot <- ggplot(data=tmp3, aes(x=(start/2+stop/2), y=dist)) + geom_line() +
                         geom_hline(yintercept=quantile(m$dist, .975, na.rm=T)) + geom_hline(yintercept=quantile(m$dist, .025, na.rm=T))

             plot_grid(freqPlot, deltaPlot, distPlot)
           }

           plotPeak(4)


  ggplot(data=m, aes(x=i, y=dist, color=chr)) + geom_line()

  o.ag.ag[J(m[grepl("7757", chr)][which.max(dist)])]

  o.ag[J(m[grepl("7757", chr)][which.max(dist)])]
