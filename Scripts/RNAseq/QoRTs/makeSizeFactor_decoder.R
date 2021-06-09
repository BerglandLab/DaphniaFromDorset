#module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0; R


#install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz",
#                   repos=NULL,
#                   type="source");

library(QoRTs)
library(data.table)
library(foreach)
library(stringr)

### set up basic decoder
  dec <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable")
  dec[,unique.ID:=gsub("pe.bam", "", id)]

  decoder <- data.table(unique.ID=gsub("pe.bam", "", dec$id),
                        lane.ID="UNKNOWN",
                        group.ID=dec$superclone,
                        sampleID=dec$rep)
  decoder[,qc.data.dir:=paste("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out", unique.ID, sep="/")]

### get the number of reads

  samps <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable")
  samps[,unique.ID:=gsub("pe.bam", "", id)]

  nreads <- foreach(samp.i=samps$SampleName, .combine="rbind")%do%{
    #samp.i <- samps[1]$unique.ID
    nReads.tot <- as.numeric(str_match(
              system(paste("grep  \"Number of input reads\" ",
                "/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/", samp.i, "_star_testLog.final.out", sep=""), intern=T),
                "[0-9]{1,}"))

    nReads.mm <- as.numeric(str_match(
              system(paste("grep  \"Number of reads mapped to multiple loci\" ",
                "/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/", samp.i, "_star_testLog.final.out", sep=""), intern=T),
                "[0-9]{1,}"))

    data.table(unique.ID=samp.i, input.read.pair.count=nReads.tot, multi.mapped.read.pair.count=nReads.mm)
  }

### merge
  decoder <- merge(decoder, nreads, by="unique.ID")


### get size factors

  res <- read.qc.results.data(infile.dir="/",
                              decoder=as.data.frame(decoder),
                              calc.DESeq2 = TRUE,
                              calc.edgeR = FALSE,
                              autodetectMissingSamples = TRUE, skip.files = c())


  sf <- get.size.factors(res)
  decoder <- merge(decoder, as.data.table(sf), by.x="unique.ID", by.y="sample.ID")
  decoder[,sample.ID:=unique.ID]

### save decoder
  write.table(decoder, sep="\t", quote=F, row.names=F,
            file="/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/decoder.txt")
