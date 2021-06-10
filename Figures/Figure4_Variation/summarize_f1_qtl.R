library(data.table)

### scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/lme4qtl_output.AxC.long.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/lme4qtl_output.CxC.long.Rdata ~/.
### scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/lme4qtl_output.v3.all.long.Rdata ~/.

### AxC
  ### load in F1 mapping data
    #load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")
    load("~/lme4qtl_output.AxC.long.Rdata")

    ### we are takign the average just because of a consequence of how the data are output. e.g., `o[id==111][term=="male"][perm==0]`
    o.ag <- o[, list(p.z=mean(p.z, na.rm=T), chisq=mean(chisq), p.aov=mean(p.aov, na.rm=T)), list(term, perm, id, chr, pos)]
    #o.ag.ag <- o.ag[,list(pr=1-mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id, chr, pos)]
    #o.ag.ag[pr==0, pr:=1/201]

    ### ANOVA output
    o.ag.perm <- o.ag[perm!=0, list(p.aov.thr=c(quantile(p.aov, .01), quantile(p.aov, .05)),
                                    q=c(.01, .05)), list(term, chr)]

    setkey(o.ag.perm, term, chr)
    setkey(o.ag, term, chr)
    o.ag.plot <- merge(o.ag[perm==0], o.ag.perm[q==0.01], allow.cartesian=T)
    o.ag.plot[,cross:="AxC"]

  ### save
    save(o.ag.plot, file="~/lme4qtl_output.AxC.long.SUMMARIZED.Rdata")


### CxC
  ### load in F1 mapping data
    #load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")
    load("~/lme4qtl_output.CxC.long.Rdata")

    ### we are takign the average just because of a consequence of how the data are output. e.g., `o[id==111][term=="male"][perm==0]`
    o.ag <- o[, list(p.z=mean(p.z, na.rm=T), chisq=mean(chisq), p.aov=mean(p.aov, na.rm=T)), list(term, perm, id, chr, pos)]
    #o.ag.ag <- o.ag[,list(pr=1-mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id, chr, pos)]
    #o.ag.ag[pr==0, pr:=1/201]

    ### ANOVA output
    o.ag.perm <- o.ag[perm!=0, list(p.aov.thr=c(quantile(p.aov, .01), quantile(p.aov, .05)),
                                    q=c(.01, .05)), list(term, chr)]

    setkey(o.ag.perm, term, chr)
    setkey(o.ag, term, chr)
    o.ag.plot <- merge(o.ag[perm==0], o.ag.perm[q==0.01], allow.cartesian=T)
    o.ag.plot[,cross:="CxC"]

  ### save
    save(o.ag.plot, file="~/lme4qtl_output.CxC.long.SUMMARIZED.Rdata")


### all
  ### load in F1 mapping data
    #load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/lme4qtl_output.AxC.long.Rdata")
    load("~/lme4qtl_output.v3.all.long.Rdata")

    ### we are takign the average just because of a consequence of how the data are output. e.g., `o[id==111][term=="male"][perm==0]`
    o.ag <- o[, list(p.z=mean(p.z, na.rm=T), chisq=mean(chisq), p.aov=mean(p.aov, na.rm=T)), list(term, perm, id, chr, pos)]
    #o.ag.ag <- o.ag[,list(pr=1-mean(p.aov[perm==0] < p.aov[perm!=0])), list(term, id, chr, pos)]
    #o.ag.ag[pr==0, pr:=1/201]

    ### ANOVA output
    o.ag.perm <- o.ag[perm!=0, list(p.aov.thr=c(quantile(p.aov, .01), quantile(p.aov, .05)),
                                    q=c(.01, .05)), list(term, chr)]

    setkey(o.ag.perm, term, chr)
    setkey(o.ag, term, chr)
    o.ag.plot <- merge(o.ag[perm==0], o.ag.perm[q==0.01], allow.cartesian=T)
    o.ag.plot[,cross:="all"]

  ### save
    save(o.ag.plot, file="~/lme4qtl_output.all.long.SUMMARIZED.Rdata")

### scp
  scp aob2x@rivanna.hpc.virginia.edu:~/lme4qtl_output.all.long.SUMMARIZED.Rdata ~/.
  scp aob2x@rivanna.hpc.virginia.edu:~/lme4qtl_output.AxC.long.SUMMARIZED.Rdata ~/.
  scp aob2x@rivanna.hpc.virginia.edu:~/lme4qtl_output.CxC.long.SUMMARIZED.Rdata ~/.
