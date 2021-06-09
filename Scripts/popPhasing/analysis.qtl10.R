#module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### makes rabbit input
  args = commandArgs(trailingOnly=TRUE)
  chr.i <- as.character(args[1])
  maxcM <- as.numeric(args[2])
  f1s.set <- as.character(args[3])
  #chr.i <- "Scaffold_1863_HRSCAF_2081"; maxcM=10; f1s.set <- "all"

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### set wd
  setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020")

### load SuperClone
  sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### which F1s?
  #f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.onlyPheno.delim")
  #f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.allF1s.delim")
  #f1s <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/F1s_to_use.all_AxC_F1s.delim")

### open GDS
  genofile <- seqOpen("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/MapJune2020_ann.seq.gds", allow.duplicate=TRUE)

### load in filter file
  snpFilter <- fread("snpsvarpulexpresentinhalf_table_20200623")


### PNPS
  ### make snp.dt
    snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                         pos=seqGetData(genofile, "position"),
                         id=seqGetData(genofile, "variant.id"),
                         numAlleles=seqNumAllele(genofile),
                         key="chr")
    setkey(snpFilter, chr, pos)
    setkey(snp.dt, chr, pos)

    snp.dt <- merge(snpFilter, snp.dt)



    seqSetFilter(genofile,
          variant.id=snp.dt$id,
          sample.id=sc[population%in%c("D8", "DBunk", "DCat", "DOil", "Dramp", "Dcat")]$clone)



          tmp <- seqGetData(genofile, "annotation/info/ANN")
          len1 <- tmp$length
          len2 <- tmp$data

          snp.dt1 <- data.table(len=rep(len1, times=len1),
                                ann=len2,
                                id=rep(snp.dt$id, times=len1))

        # Extracting data between the 2nd and third | symbol
          snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
          snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]


        # Collapsing additional annotations to original SNP vector length
          snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), col2=paste(gene, collapse=",")),
                                list(variant.id=id)]

          snp.dt1.an[,class:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
          snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$col2,"\\,")[[1]]]

          pnps <- snp.dt1.an[,list(NS=sum(class=="missense_variant"), Syn=sum(class=="synonymous_variant")), list(gene)]


          save(pnps, file="~/pnps.Rdata")

#### DNDS
  load("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/dpfiltsnps_20200623.Rdata")

  seqSetFilter(genofile,
        variant.id=dpfiltsnps$variant.ids,
        sample.id="2018_Pulicaria_Pond21_22")

                  tmp <- seqGetData(genofile, "annotation/info/ANN")
                  len1 <- tmp$length
                  len2 <- tmp$data

                  snp.dt1 <- data.table(len=rep(len1, times=len1),
                                        ann=len2,
                                        id=rep(dpfiltsnps$variant.ids, times=len1))

                # Extracting data between the 2nd and third | symbol
                  snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
                  snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]


                # Collapsing additional annotations to original SNP vector length
                  snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), col2=paste(gene, collapse=",")),
                                        list(variant.id=id)]

                  snp.dt1.an[,class:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
                  snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$col2,"\\,")[[1]]]

                  dnds <- snp.dt1.an[,list(NS=sum(class=="missense_variant"), Syn=sum(class=="synonymous_variant")), list(gene)]

          save(dnds, pnps, file="~/dnds_pnps.Rdata")



scp aob2x@rivanna.hpc.virginia.edu:~/dnds_pnps.Rdata ~/.
  load(file="~/dnds_pnps.Rdata")

  setnames(dnds, c("NS", "Syn"), c("dn", "ds"))
  setnames(pnps, c("NS", "Syn"), c("pn", "ps"))

          m <- merge(dnds, pnps, by="gene")


          m[,NI:=log((pn/ps)/(dn/ds))]

          pnps[,r:=NS/Syn]
          pnps[,n:=NS+Syn]
          m[gene=="Daphnia00787"]

          m.ag <- m[dn+ds>0 & pn+ps>0,list(p=chisq.test(matrix(c(pn, ps, dn, ds), nrow=2, byrow=T))$p.value), list(gene)]

hist(m[n>5][NI!=Inf]$NI)
