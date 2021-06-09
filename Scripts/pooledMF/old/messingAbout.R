### libraries
  library(data.table)
  library(ggplot2)
  library(doMC)
  registerDoMC(20)
  library(cowplot)
  library(poolfstat)    #; install.packages("poolfstat")

### data
  load("/mnt/sammas_storage/bergland-lab/alan/totrdfilt.Rdata")
  load("/mnt/sammas_storage/bergland-lab/alan/totalADRDlongall.Rdata")

### generate a SNP metatable with filtering info in there
  snp.dt <- unique(totrdfilt[,c("chr", "pos"),with=F])
  snp.dt[,karen_rd_pass:=T]

  snp.dt.ag <- snp.dt[,.N, chr]
  snp.dt.ag[,chr_pass:=N>15000]

  setkey(snp.dt.ag, chr)

### merge in filtering into long dataset
  setkey(geno, chr, pos)
  setkey(snp.dt, chr, pos)

  geno <- merge(geno, snp.dt, all.x=T)

  setkey(geno, chr)
  setkey(snp.dt.ag, chr)

  geno <- merge(geno, snp.dt.ag, all.x=T)
  geno[is.na(karen_rd_pass), karen_rd_pass:=FALSE]
  geno[is.na(chr_pass), chr_pass:=FALSE]

  table(geno$karen_rd_pass); table(is.na(geno$karen_rd_pass))
  table(geno$chr_pass); table(is.na(geno$chr_pass))

  pwFst <- function(m, f, m.n, f.n) {
    # m <- "D8Male1"; f <- "D8PE1"; m.n=35; f.n=50

    ### make into poolfstat object
      res <- new("pooldata")
      res@npools = 2
      res@nsnp = length(geno[ref!="N"][karen_rd_pass==T][Sample=="D8Male1"]$RD)

      res@refallele.readcount = cbind(geno[ref!="N"][karen_rd_pass==T][Sample==m]$RD,
                                      geno[ref!="N"][karen_rd_pass==T][Sample==f]$RD)

      res@readcoverage = cbind(geno[ref!="N"][karen_rd_pass==T][Sample==m]$TD,
                               geno[ref!="N"][karen_rd_pass==T][Sample==f]$TD)

      res@snp.info = cbind(geno[ref!="N"][karen_rd_pass==T][Sample=="D8Male1"]$chr,
                           geno[ref!="N"][karen_rd_pass==T][Sample=="D8Male1"]$pos,
                           geno[ref!="N"][karen_rd_pass==T][Sample=="D8Male1"]$ref,
                           geno[ref!="N"][karen_rd_pass==T][Sample=="D8Male1"]$alt)

      res@poolsizes = 2*c(m.n, f.n)
      res@poolnames = c(m, f)

    ## Estimate Fst
      o <- computeFST(res, method="Anova")

      o <- data.table(chr=geno[ref!="N"][karen_rd_pass==T][Sample=="D8Male1"]$chr,
                      pos=geno[ref!="N"][karen_rd_pass==T][Sample=="D8Male1"]$pos,
                      fst=o$snp.FST,
                      chr_pass=geno[ref!="N"][Sample=="D8Male1"][karen_rd_pass==T]$chr_pass)


      o$chr_use <- o$chr
      o[chr_pass==F, chr_use:="extra_stuff"]
      o <- o[!is.na(fst) & !is.nan(fst)]
      o[, r := rank(fst)]

      o[, z:= qnorm(r/(length(fst)+1), 0, 1)]

    ### return
      o[,m:=m]
      o[,f:=f]

      return(o)
  }

### some basic contrasts
  fst.o <- rbind(pwFst(m="D8Male1", f="D8PE1", m.n=35, f.n=50),
                 pwFst(m="D8Male2", f="D8PE2", m.n=35, f.n=50),
                 pwFst(m="D8Male1", f="D8PE2", m.n=35, f.n=50),
                 pwFst(m="D8Male2", f="D8PE1", m.n=35, f.n=50),
                 pwFst(m="D8Male1", f="D8Male2", m.n=35, f.n=35),
                 pwFst(m="D8PE1", f="D8PE2", m.n=50, f.n=50),
                 pwFst(m="DBunkMale", f="DBunkPE1", m.n=100, f.n=50),
                 pwFst(m="DBunkMale", f="DBunkPE2", m.n=100, f.n=50),
                 pwFst(m="DBunkPE1", f="DBunkPE2", m.n=50, f.n=50),
                 pwFst(m="DCatMale", f="DCatPE1", m.n=100, f.n=50),
                 pwFst(m="DCatMale", f="DCatPE2", m.n=100, f.n=50),
                 pwFst(m="DCatPE1", f="DCatPE2", m.n=50, f.n=50))


### basic plots
  fst.o[,pair:=paste(m, f, sep="_")]
  ggplot(data=fst.o, aes(x=pair, y=fst)) + geom_violin() + coord_flip()


### sliding window

  ### make windows

    d.ag <- geno[,list(n=length(pos), start=min(pos), stop=max(pos)), list(chr)]
    step.size <- 5000
    window.size <- 50*step.size

    wins <- foreach(chr.i=d.ag$chr, .combine="rbind")%dopar%{
      print(chr.i)
      if( (d.ag[chr==chr.i]$stop - d.ag[chr==chr.i]$start)>window.size) {
        tmp <- data.table(chr=chr.i, start=seq(from=d.ag[chr==chr.i]$start, to=d.ag[chr==chr.i]$stop - window.size, by=step.size))
        tmp[,stop:=start+window.size]
      } else {
        tmp <- data.table(chr=chr.i, start=d.ag[chr==chr.i]$start, stop=d.ag[chr==chr.i]$stop)
      }
      tmp
    }
    wins[,i:=1:dim(wins)[1]]

  ### iterate through window

    setkey(fst.o, chr, pos)
    fst.o.ag <- foreach(i=c(1:dim(wins)[1]), .combine="rbind", .errorhandling="remove")%dopar%{
      print(i)
      ### get data; i<-10887
        tmp <- fst.o[J(data.table(chr=wins$chr[i], pos=wins$start[i]:wins$stop[i], key="chr,pos")), nomatch=0]

        #ggplot(data=tmp[ref!="N"][karen_rd_pass==T], aes(x=pos, y=propalt, group=Sample, color=sex)) + geom_line() +facet_grid(pond~.)

      ### might be clunky at first
        tmp[,list(i=i,
                  mu.z=mean(z, na.rm=T), med.z=median(z, na.rm=T),
                  mu.fst=mean(fst, na.rm=T), med.fst=median(fst, na.rm=T),
                  n.snps=length(z[!is.na(z)]),
                  p.z=t.test(z)$p.value),
             list(m, f)]

    }



    fst.o.ag[,pair:=paste(m, f, sep="_")]
    setkey(fst.o.ag, i)
    setkey(wins, i)
    fst.o.ag <- merge(fst.o.ag, wins)

    setkey(fst.o.ag, chr)
    setkey(d.ag, chr)
    fst.o.ag <- merge(fst.o.ag, d.ag)
    fst.o.ag[,pond:="foo"]
    fst.o.ag[grepl("D8", pair),pond:="D8"]
    fst.o.ag[grepl("DBunk", pair),pond:="DBunk"]
    fst.o.ag[grepl("DCat", pair),pond:="DCat"]

    save(fst.o.ag, fst.o, wins, d.ag, geno, file="~/male_female_fst.Rdata")

### plot
    load(file="~/male_female_fst.Rdata")
    fst.o.ag[,pond:="foo"]
    fst.o.ag[grepl("D8", pair),pond:="D8"]
    fst.o.ag[grepl("DBunk", pair),pond:="DBunk"]
    fst.o.ag[grepl("DCat", pair),pond:="DCat"]

    pond.i <- "D8"
    ggplot(data=fst.o.ag[n.x>100][pond==pond.i], aes(x=i, y=mu.fst, color=chr)) +
    geom_vline(xintercept=fst.o.ag[n.x>100][pair=="D8Male1_D8PE2"][which.max(med.fst)][pond==pond.i]$i) +
    geom_vline(xintercept=fst.o.ag[n.x>100][pair=="D8Male2_D8PE2"][which.max(med.fst)][pond==pond.i]$i) +
    geom_line() +
    geom_point(data=fst.o.ag[n.x>100][med.fst>max(fst.o.ag[pair=="D8PE1_D8PE2"][n.x>100]$med.fst)][pond==pond.i], aes(x=i, y=mu.fst, color=chr)) +
    facet_wrap(~pair, ncol=1) +
    theme(legend.position = "none")


    fst.o.ag[n.snps>100][which.max(med.fst)]




    ggplot(data=fst.o.ag[n.x>100][pond==pond.i][chr==fst.o.ag[n.x>100][pond==pond.i][which.max(med.fst)]$chr],
            aes(x=(start.x/2 + stop.x/2), y=med.fst, color=chr)) +
    geom_vline(xintercept=fst.o.ag[n.x>100][pond==pond.i][pair=="D8Male1_D8PE2"][which.max(mu.fst)]$i) +
    geom_vline(xintercept=fst.o.ag[n.x>100][pond==pond.i][pair=="D8Male2_D8PE2"][which.max(mu.fst)]$i) +
    geom_line() +
    geom_point(data=fst.o.ag[n.x>100][n.y>0][med.fst>max(fst.o.ag[pair=="D8PE1_D8PE2"][n.x>100]$med.fst)][chr==fst.o.ag[n.x>100][n.y>0][which.max(med.fst)]$chr],
          aes(x=(start.x/2 + stop.x/2), y=med.fst, color=chr)) +
    facet_wrap(~pair, ncol=1) +
    theme(legend.position = "none")







    setkey(geno, chr, pos)
    #setkey(wins, chr, pos)
    i <- fst.o.ag[n.x>100][pond=="D8"][pair=="D8Male2_D8PE2"][which.max(mu.fst)]$i
    buffer <- 0
    tmp <- geno[J(data.table(chr=wins$chr[i], pos=(wins$start[i]-buffer):(wins$stop[i]+buffer), key="chr,pos")), nomatch=0]

    ggplot() +
    geom_boxplot(data=tmp, aes(y=propalt, x=Sample, fill=sex, group=Sample)) +
    facet_grid(pond~., scales="free_y") + coord_flip()

    ggplot() +
    geom_density(data=tmp, aes(propalt, color=sex, group=Sample)) +
    facet_grid(pond~., scales="free_y")


    fst.o.ag[pair=="D8Male1_D8PE2"][n.x>100][med.fst>.105]




  fst.o.ag[n>100][m=="D8Male1"]


  plot(-log10(p)~i, o.ag[n>100])


  totalADRDlongall.Rdata
