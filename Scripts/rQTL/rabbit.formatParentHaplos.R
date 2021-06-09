#module load intel/18.0 intelmpi/18.0 R/3.6.3; R


### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(SeqArray)

### make SNP table
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       pos=seqGetData(genofile, "position"),
                       id=seqGetData(genofile, "variant.id"),
                       numAlleles=seqNumAllele(genofile))
  setkey(snp.dt, id)


  loadDat <- function(fn) {
    print(fn)
    #fn <- fns[12]
    pp <- fread(fn, skip=1, nrows=6, header=T, fill=T)
    ppl <- melt(pp, id.vars="marker")

    setnames(ppl, c("marker", "variable"), c("allele", "id"))

    setkey(ppl, "allele")

    ppl <- ppl[J(c("A_Maternal", "A_Paternal", "C_Maternal", "C_Paternal"))]
    ppl[,id:=as.numeric(as.character(id))]

    setkey(ppl, "id")
    setkey(snp.dt, "id")
    ppl <- merge(ppl, snp.dt)
    return(ppl)
  }

  fns <- system("ls /scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/*/*out.post.csv", intern=T)
  parental <- foreach(x=fns)%do%loadDat(x)
  parental <- rbindlist(parental)

  save(parental, file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/parental.Rdata")
