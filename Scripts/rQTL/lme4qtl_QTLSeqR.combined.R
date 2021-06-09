### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(patchwork)


### Pool-Seq

  load("~/gprime_peaks.replicates.Rdata")
  load("~/gprime_peaks.groups.250K.05.Rdata")

  setnames(peaks, "CHROM", "chr")
  setnames(gprime, "CHROM", "chr")
  setnames(peaks.groups, "CHROM", "chr")
  setnames(gprime.groups, "CHROM", "chr")

  peaks <- rbind(peaks, peaks.groups, fill=T)
  gprime <- rbind(gprime, gprime.groups, fill=T)

  peaks[,rep:=gsub("NA", "", paste(rep, group, sep=""))]
  gprime[,rep:=gsub("NA", "", paste(rep, group, sep=""))]

  poolseq <- ggplot() +
  geom_vline(data=peaks[is.na(group)], aes(xintercept=posMaxGprime), color="black") +
  geom_line(data=gprime[is.na(group)], aes(x=POS, y=Gprime, color=chr)) +
  facet_grid(rep~chr, scales="free_x") +
  theme(legend.position = "none")


### F1
  load("~/lme4qtl_output.AxC.long.Rdata")

  o.ag <- o[, list(p.z=mean(p.z, na.rm=T), chisq=max(chisq), p.aov=min(p.aov, na.rm=T)), list(term, perm, id, chr, pos)]
  o.ag.ag <- o.ag[,list(pr=mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id, chr, pos)]
  o.ag.ag[pr==1, pr:=1/201]

  ### ANOVA output
  o.ag.perm <- o.ag[perm!=0, list(p.aov=c(quantile(p.aov, .01), quantile(p.aov, .05)),
                                  q=c(.01, .05)), list(term, chr)]



  f1.plot <- ggplot() +
  geom_vline(data=peaks[is.na(group)], aes(xintercept=posMaxGprime, linetype=as.factor(rep))) +
  geom_line(data=o.ag[perm==0], aes(x=pos, y=-log10(p.aov), color=chr)) +
  geom_hline(data=o.ag.perm, aes(yintercept=-log10(p.aov), linetype=as.factor(q))) +
  facet_grid(term~chr, scales="free_x") +
  theme(legend.position = "none")

### combined plot
  combined <- poolseq / f1.plot

  ggsave(combined, file="~/combined_mapping.png")

### overlap test
  setkey(o.ag, perm )



### pos-max
  peaks[,pos:=posMaxGprime]
  setkey(peaks, chr, pos)

  setkey(o.ag, chr, pos)
  m <- merge(peaks, o.ag, all.x=T)
  m[perm==0]






#### this stretches peaks
  peaks.long <- foreach(peaks.i=1:dim(peaks)[1], .combine="rbind")%do%{
    data.table(peaks.i=peaks.i, rep=peaks$rep[peaks.i], chr=peaks$chr[peaks.i], pos=peaks$start[peaks.i]:peaks$end[peaks.i])
  }
  setkey(peaks.long, chr, pos)








  overlap <- foreach(perm.i=c(0:200), .combine="rbind", .errorhandling="remove")%do%{
    #perm.i<-0
    message(perm.i)
    tmp <- o.ag[J(perm.i)]
    setkey(tmp, chr, pos)

    tmp <- merge(tmp, peaks.long, all.x=T)

    foreach(term.i=c("male", "fill"), .combine="rbind", .errorhandling="remove")%do%{
      foreach(rep.i=c(1,2), .combine="rbind", .errorhandling="remove")%do%{
        mat <- table(tmp[term==term.i]$p.aov <= mean(o.ag.perm[term==term.i][q==.01]$p.aov),
                     !is.na(tmp[term==term.i]$peaks.i))
        ft <- fisher.test(mat)

        data.table(or=ft$estimate, p=ft$p.value, perm=perm.i, term=term.i)
      }
    }

  }

  setkey(o.ag, perm)
  maxstat <- foreach(perm.i=c(0:200), .combine="rbind", .errorhandling="remove")%do%{
    #perm.i<-0
    message(perm.i)
    tmp <- o.ag[J(perm.i)]
    setkey(tmp, chr, pos)

    tmp <- merge(tmp, peaks.long, all.x=T)

    tmp[,list(min.p=min(p.aov)), list(rep, peaks.i, term, perm)]


  }

  maxstat

  ggplot(data=maxstat[sort(perm, decreasing=T)], aes(x=peaks.i, y=-log10(min.p), color=as.factor(perm==0))) +
  geom_point() +
  facet_grid(term~rep)


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
