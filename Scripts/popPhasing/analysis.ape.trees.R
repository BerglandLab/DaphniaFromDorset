#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

#scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/popPhase/trees/region.fasta ~/.

### libraries
  library(ape)
  library(data.table)

### get file name from input stream

  args = commandArgs(trailingOnly=TRUE)
  fasta.fn=args[[1]]
  out.fn=args[[2]]

  #fasta.fn <- "/scratch/aob2x/daphnia_hwe_sims/popPhase/trees/region.fasta"
  #fasta.fn <- "/dev/shm/aob2x/1/1/Scaffold_2217_HRSCAF_2652_5173221_5223221.fasta"

### load fasta as DNAbin
  message("load ing fasta")
  db <- read.FASTA(fasta.fn, type = "DNA")

### dist
  message("get distances")
  gd <- dist.dna(db, as.matrix=T)
  gdm <- as.matrix(gd)
  dim(gdm)
  gdml <- as.data.table(expand.grid(gdm)[,1])
  gdml

### tree
  message("making tree")

  njo <- bionj(gd)
  #yBoots <- boot.phylo(njo, as.matrix(db), function(xx) bionj(dist.dna(xx)), B = 10,  mc.cores = 1)


  #njo <- root(njo, outgroup=njo$tip.label[grepl("March20_2018_DBunk_38.1.fa", njo$tip.label)][1])

### make info object
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


  save(cdl, gd, njo, d, file=out.fn)
