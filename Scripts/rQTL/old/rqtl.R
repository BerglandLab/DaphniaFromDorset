### libraries
  library(qtl)
  library(data.table)
  library(ggplot2)
  library(patchwork)

############################
#### Prep the input data ###
############################

### load data and convert: [M]arkers, [H]eader
  f1s <- fread("~/AxC_F1.csv") ### Comes from `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbitRedo/plotRabbit.R`
  f1s[,marker:=paste(chr.x, id, sep="_")]
  f1s[,diplo:=as.numeric(as.factor(founder))]
  setkey(f1s, marker)

  tmp <- data.table(marker=sample(unique(f1s$marker), 6000, replace=F))
  setkey(tmp, marker)
  f1s.sub <- f1s[J(tmp)]
  setkey(f1s.sub, id)

  markers<- dcast(f1s.sub , clone ~ id, value.var=list("diplo"))
  markers[1:5,1:4]


  chrs <- f1s.sub[,list(chr=unique(chr.x), pos=unique(pos)), id]

  header <- as.data.table(rbind(c("", chrs$chr), c("", chrs$pos)))
  setnames(header, names(header), names(markers))

  markers[1:5,1:4]
  header[1:2,1:4]

  mh <- rbind(header, markers)
  mh[1:5,1:4]

### set up [P]henotype data
  ### This uses the raw data
    # load(file="~/m3epp.Rdata")
    # phenos <- m3epp.ag[Type=="AxCF1", c("mu.epp", "mu.fill", "Clone"), with=F]

    # phenos[Clone=="D818111", Clone:="April5_2018_D8_18111"]
    # phenos[Clone=="D818106", Clone:="April5_2018_D8_18106"]
    # phenos[Clone=="D818025", Clone:="March20_2018_D8_18025"]
    # phenos[Clone=="D818028", Clone:="March20_2018_D8_18028"]
    # phenos[Clone=="D818030", Clone:="March20_2018_D8_18030"]
    # phenos[Clone=="D818010", Clone:="March20_2018_D8_18010"]

    # setnames(phenos, "Clone", "clone")

  ### This uses the BLUP phenotypes, extracted from the rQTL file that Karen had made
    AxCF1 <- read.cross("csv","","~/AxCF1genoandphenoreconstruct_sub3.csv", genotypes=NULL)
    AxCF1$pheno

    phenos <- as.data.table(AxCF1$pheno)
    setnames(phenos, "CloneID", "clone")

    ### ---> Tack in new phenotype data here

### merge to make rQTL file
  setkey(phenos, "clone")
  setkey(mh, "clone")

  mhp <- merge(phenos, mh, all.y=T)

  setcolorder(mhp, c("eppresid", "embresid", "SCB", "clone", names(markers)[-1]))

  mhp[1:5,1:4]
  dim(mhp)

  write.csv(mhp, file="~/mhp.csv", quote=F, row.names=F, na="")

### Run rQTL
  ### read cross object
    AxCF1 <- read.cross("csv","","~/mhp.csv", crosstype="4way", genotypes=NULL)

  ### marker regression

    mr1 <- scanone(AxCF1, pheno.col=1, method="mr")
    mr2 <- scanone(AxCF1, pheno.col=2, method="mr")

    mr1.dt <- as.data.table(mr1)
    mr2.dt <- as.data.table(mr2)

### Load pooled data
    load("~/peaks.Rdata")
    peaks <- fread("/Users/alanbergland/peaks.csv")

    setnames(peaks, "CHROM", "chr")
    setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))

     Scaffold_7757_HRSCAF_8726 8706276

### plot
    eppresid.plot <- ggplot() +
    #geom_hline(yintercept=summary(perm)[2]) +
    geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
    geom_line(data=mr1.dt, aes(x=pos, y=lod, color=chr), size=1) +
    facet_grid(.~chr, scales="free_x") +
    theme(strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
    ggtitle("epp.resid")


    #setnames(gprime, c("CHROM", "POS"), c("chr", "pos"))
    Gprime.plot <- ggplot(data=gprime, aes(x=pos, y=Gprime, color=chr)) +
    geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
    geom_line(size=.75) +
    facet_grid(.~chr, scales="free_x") +
    theme(strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
    ggtitle("pooledWild")

    embresid.plot <- ggplot() +
    #geom_hline(yintercept=summary(perm)[2]) +
    geom_vline(data=peaks, aes(xintercept=posMaxGprime), linetype="dashed") +
    geom_line(data=mr2.dt, aes(x=pos, y=lod, color=chr), size=1) +
    facet_grid(.~chr, scales="free_x") +
    theme(strip.text.x = element_text(size = 8),
          axis.text.x = element_text(angle = 90, size=.5), legend.position = "none") +
    ggtitle("emb.resid")


### final plot
  eppresid.plot / Gprime.plot / embresid.plot


### save
  tar czvf rQTL.inputFiles.tar.gz \
  ~/AxC_F1.csv \
  ~/AxCF1genoandphenoreconstruct_sub3.csv \
  ~/peaks.Rdata \
  /Users/alanbergland/peaks.csv
