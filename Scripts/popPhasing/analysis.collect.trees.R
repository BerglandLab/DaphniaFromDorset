#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


### libraries
  library(data.table)
  library(foreach)
  library(ape)
  library(doMC)
  registerDoMC(20)

### input files: genome-wide distribution


  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/popPhase/trees250K/", "Rdata", full.names=T)
  length(fn)

  cdl.list <- foreach(i=fn)%dopar%{
    message(i)
    #i <- fn[1]
    load(i)
    #return(cdl)
    cdl[,window:=tstrsplit(i1, ";")[[2]]]

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
    #cdl[,cd_bin:=factor(cd_bin, seq(from=0, to=.1, by=.001))]

    njo$tip.label <- gsub(".fa", "", tstrsplit(njo$tip.label, ";")[[1]])

    njo <- root(njo, njo$tip.label[grepl("obtusa", njo$tip.label)])

    o.tmp <- list(cdl[,list(n=.N), list(sp.group, pond.group, window, cd_bin)],
                  njo)

    names(o.tmp) <- c(cdl[1]$window, cdl[1]$window)

    o.tmp
  }

  cdl.o <- lapply(cdl.list, function(x) x[[1]])
  cdl.o <- rbindlist(cdl.o)
  cdl.genome <- cdl.o[,list(n=sum(n)), list(sp.group, pond.group, cd_bin)][!is.na(sp.group) & !is.na(pond.group)]


  cdl.tree <- lapply(cdl.list, function(x) x[[2]])
  names(cdl.tree) <- lapply(cdl.list, function(x) names(x)[2])


### qtl
  wins <- fread("/scratch/aob2x/daphnia_hwe_sims/popPhase/windows.delim")

  target.fn <- sapply(tail(wins, 14)$V2, function(x) fn[grepl(x, fn)])

  length(target.fn)
  cdl.qtl <- foreach(i=target.fn)%dopar%{
    message(i)
    #i <- fn[1]
    load(i)
    #return(cdl)
    cdl[,window:=tstrsplit(i1, ";")[[2]]]

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


    #cdl[,list(cd_mean=mean(cd), cd_sd=sd(cd)), list(sp.group, pond.group, window)]
  }
  cdl.qtl <- rbindlist(cdl.qtl)


### save


  save(cdl.genome, cdl.qtl, cdl.o, cdl.tree, cdl.list, file="~/cdlo_250K.Rdata")
