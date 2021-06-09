#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

#scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/popPhase/trees/region.fasta ~/.

### libraries
  library(ape)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(5)
  
### get file name from input stream

  args = commandArgs(trailingOnly=TRUE)
  fasta.fn=args[[1]]
  out.fn=args[[2]]

  #fasta.fn <- "/scratch/aob2x/daphnia_hwe_sims/popPhase/trees/region.fasta"
  #fasta.fn <- "/dev/shm/aob2x/1/1/Scaffold_2217_HRSCAF_2652_5173221_5223221.fasta"

### load fasta as DNAbin
  message("load ing trees")
  #
  bs <- read.tree(paste(fasta.fn, ".boottrees", sep=""))

### make info object
  cdl <- foreach(i=c(1:length(bs)))%dopar%{
    message(i)
    njo <- bs[[i]]

    d <- data.table(tip.label=njo$tip.label)
    d[,pond:=tstrsplit(tip.label, "_")[[3]]]
    d[,sample.id:=tstrsplit(tip.label, "\\.")[[1]]]

    #d[,A:=ifelse(grepl("April_2017_D8_213", label), "A", "")]
    #d[,C:=ifelse(grepl("April_2017_D8_151", label), "C", "")]
    d[pond%in%c("D8", "DBunk", "DCat", "DOil", "Dramp", "Dcat"),group:="DWT"]
    d[pond=="D10", group:="D10"]
    d[grepl("W", pond), group:="W"]
    d[grepl("Pond22", pond), group:="pulicaria"]
    d[grepl("Pond22", pond), group:="pulicaria"]
    d[grepl("March20_2018_DBunk_38", tip.label), group:="obtusa"]
    d[grepl("April_2017_D8_213", tip.label), group:="A"]
    d[grepl("April_2017_D8_151", tip.label), group:="C"]
    d[group%in%c("pulicaria", "obtusa"), species:=group]
    d[!group%in%c("pulicaria", "obtusa"), species:="pulex"]


    ### cd
    cd <- cophenetic(njo)
    cd[diag(cd)] <- NA
    cd[upper.tri(cd)] <- NA
    cdl <-data.table(cd=expand.grid(cd)[,1],
                      i1=rep(rownames(cd), length(rownames(cd))),
                      i2=rep(colnames(cd), each=length(rownames(cd))))
    cdl <- na.omit(cdl)

    cdl[,clone1:=tstrsplit(i1, "\\.")[[1]]]
    cdl[,clone2:=tstrsplit(i2, "\\.")[[1]]]

    setnames(cdl, "i1", "tip.label")
    cdl <- merge(cdl, d, by="tip.label",)
    setnames(cdl, "tip.label", "i1")

    setnames(cdl, "i2", "tip.label")
    cdl <- merge(cdl, d, by="tip.label")
    setnames(cdl, "tip.label", "i2")


    cdl[species.x=="pulex" & species.y=="pulex", sp.group:="pulex-pulex"]
    cdl[(species.x=="pulex" & species.y=="pulicaria") | (species.y=="pulex" & species.x=="pulicaria"), sp.group:="pulex-pulicaria"]
    cdl[(species.x=="pulex" & species.y=="obtusa") | (species.y=="pulex" & species.x=="obtusa"), sp.group:="pulex-obtusa"]
    cdl[,sp.group:=factor(sp.group, levels=c("pulex-pulex", "pulex-pulicaria", "pulex-obtusa"))]


    cdl[(group.x=="DWT" & group.y=="DWT"), pond.group:="DWT-DWT"]
    cdl[(group.x=="DWT" & group.y=="D10") | (group.y=="DWT" & group.x=="D10"), pond.group:="DWT-D10"]
    cdl[(group.x=="DWT" & group.y=="W") | (group.y=="DWT" & group.x=="W"), pond.group:="DWT-W"]
    cdl[sp.group!="pulex-pulex", pond.group:="all"]
    cdl[group.x%in%c("A", "C") & group.y%in%c("A", "C"), pond.group:="DWT-DWT"]
    cdl[,pond.group:=factor(pond.group, levels=c("DWT-DWT", "DWT-D10", "DWT-W", "all"))]



    cdl[,window:=last(tstrsplit(i1, "_"))]

    ### fix bug
      cdl[sample.id.x=="pulicaria", group.x:=   "pulicaria"]
      cdl[sample.id.x=="pulicaria", species.x:= "pulicaria"]
      cdl[sample.id.x=="pulicaria", pond.x:=    "pulicaria"]

      cdl[sample.id.y=="pulicaria", group.y:=   "pulicaria"]
      cdl[sample.id.y=="pulicaria", species.y:= "pulicaria"]
      cdl[sample.id.y=="pulicaria", pond.y:=    "pulicaria"]

      cdl[sample.id.x=="obtusa", group.x:=   "pulicaria"]
      cdl[sample.id.x=="obtusa", species.x:= "pulicaria"]
      cdl[sample.id.x=="obtusa", pond.x:=    "pulicaria"]

      cdl[sample.id.y=="obtusa", group.y:=   "obtusa"]
      cdl[sample.id.y=="obtusa", species.y:= "obtusa"]
      cdl[sample.id.y=="obtusa", pond.y:=    "obtusa"]

      cdl[(species.x=="pulex" & species.y=="pulicaria") | (species.y=="pulex" & species.x=="pulicaria"), sp.group:="pulex-pulicaria"]
      cdl[(species.x=="pulex" & species.y=="obtusa") | (species.y=="pulex" & species.x=="obtusa"), sp.group:="pulex-obtusa"]

      cdl[(species.x=="pulicaria" & species.y=="obtusa") | (species.y=="pulicaria" & species.x=="obtusa"), sp.group:="pulicaria-obtusa"]
      cdl[sp.group!="pulex-pulex", pond.group:="all"]

      cdl[group.x%in%c("A", "C") & group.y%in%c("A", "C"), pond.group:="A-C"]

      cdl <- cdl[sample.id.x!=sample.id.y]



      cdl[,cd_bin:=round(cd, 5)]

      cdl[is.na(window), window:=unique(na.omit(cdl$window))]
      cdl.tmp <- cdl[,list(n=.N), list(sp.group, pond.group, window, cd_bin)]
      cdl.tmp[,boot:=i]
      cdl.tmp
    }
  cdl <- rbindlist(cdl)

  njo <- read.tree(paste(fasta.fn, ".bionj", sep=""))
  bs <- read.tree(paste(fasta.fn, ".boottrees", sep=""))




  save(cdl, njo, bs, file=out.fn)
