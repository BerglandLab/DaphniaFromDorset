#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### library
  library(data.table)
  library(foreach)
  library(SeqArray)

### load AxC F1 reconstruction (from `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.convert_output.R`)
  AxC <- fread(file="/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all_AxC.csv")

  snps.dt <- AxC[,list(id=unique(id)), list(chr=chr.y, pos=pos)]

### load CxC genotypes at these positions

  ### set wd
    setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

  ### load SuperClone
    sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

  ### open GDS
    genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

  ### pull out genotypes for CxC @ snp set
    seqSetFilter(genofile,
                  variant.id=snps.dt$id,
                  sample.id=c(sc[SC=="C"][which.max(medrd)]$clone,
                              sc[OneLiterPheno==1][AxCF1Hybrid==0][SC=="selfedC"]$clone))

    genomat <- as.data.table(t(seqGetData(genofile, "$dosage")))
    setnames(genomat, seqGetData(genofile, "sample.id"))
    setnames(genomat, sc[SC=="C"][which.max(medrd)]$clone, "C")
    genomat[,id:=seqGetData(genofile, "variant.id")]

    ### are CxC really looking as they should? yes
      table(expand.grid(as.matrix(genomat[C==0]))[,1])
      table(expand.grid(as.matrix(genomat[C==1]))[,1])
      table(expand.grid(as.matrix(genomat[C==2]))[,1])


































  ### make long
    gml <- melt(genomat, id.vars="id")
    setnames(gml, "variable", "clone")
    gml[value==0,founder:="E"]
    gml[value==1,founder:="F"]
    gml[value==2,founder:="G"]


  ### merge w/ AxC
    allF1 <- rbind(AxC, gml, fill=T)

    length(unique(allF1$clone))

### load in phenotype data
    load(file="~/F1_pheno.Rdata")
    mm[is.na(clone), clone:=SCB]
    setkey(mm, clone)

### iterate through
  setkey(allF1, id)
  o <- foreach(i=snps.dt$id[1:1000])%do%{
    if(which(i==snps.dt$id)%%100==0) print(paste(which(i==snps.dt$id)))

    #i=snps.dt[pos==6229430 & chr=="Scaffold_9199_HRSCAF_10755"]$id
    #i<-15618
    #i<-3297
    tmp.geno <- allF1[J(i)]
    tmp.gp <- merge(tmp.geno, mm, by="clone")


    #tmp.lm <- lm(fill~founder, tmp.gp)
    #summary(tmp.lm)
    #
    tmp.aov <- summary(aov(fill~founder, tmp.gp))

    data.table(id=i,
               F=tmp.aov[[1]]$F[1],
               p=tmp.aov[[1]]$Pr[1])
  }
  o <- rbindlist(o)
