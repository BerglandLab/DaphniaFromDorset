#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### library
  library(data.table)
  library(foreach)
  library(stringr)
  library(GenomicFeatures)
  library(patchwork)

### library size data

  samps <- fread("/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/sampleTable")
  samps[,unique.ID:=gsub("pe.bam", "", id)]

  nreads <- foreach(samp.i=samps$SampleName, .combine="rbind")%do%{
    #samp.i <- samps[1]$unique.ID
    nReads <- as.numeric(str_match(
              system(paste("grep  \"Number of input reads\" ",
                "/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/", samp.i, "_star_testLog.final.out", sep=""), intern=T),
                "[0-9]{1,}"))
    data.table(samp=samp.i, nReads=nReads)
  }

### load data
  fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out", full.names=T)

  genebody <- foreach(fn.i=fn, .errorhandling="remove")%do%{
    #fn.i <- fn[1]
    tmp <- fread(paste(fn.i, "QC.geneBodyCoverage.genewise.txt.gz", sep="/"), header=T)

    gbl <- melt(tmp, "GENE_ID")
    gbl.ag <- gbl[,list(aveExp=mean(value)), list(GENE_ID)]
    gbl <- merge(gbl, gbl.ag, by="GENE_ID")
    gbl[,samp:=last(tstrsplit(fn.i, "/"))]
    gbl <- merge(gbl, nreads, by="samp")

  }
  genebody <- rbindlist(genebody)
  genebody <- merge(genebody, samps, by.x="samp", by.y="SampleName")

### get some features of these genes

  save(genebody, file="~/genebody.Rdata")


#### on MBP
# scp aob2x@rivanna.hpc.virginia.edu:~/genebody.Rdata ~/.

### libraries
  library(data.table)
  library(ggplot2)
  library(GenomicFeatures)


### get gene features
  txdb <- makeTxDbFromGFF(file="/Users/alanbergland/daphnia_ref/Daphnia.aed.0.6.gff", format="gff3")
  tx.len <- transcriptLengths(txdb, with.cds_len=TRUE,
                                    with.utr5_len=TRUE,
                                    with.utr3_len=TRUE)
  tx.len <- as.data.table(tx.len)
  tx.len.ag <- tx.len[,list(nexon=max(nexon), cds_len=max(cds_len)), list(GENE_ID=gene_id)]

### load data
  load("~/genebody.Rdata")

  genebody[,exp.norm:=(value/nReads)/(aveExp/nReads)]
  genebody[,var:=as.numeric(as.character(variable))]

### merge
  genebody <- merge(genebody, tx.len.ag, by="GENE_ID")
  genebody.ag <- genebody[,list(mu=mean(exp.norm, na.rm=T), sd=sd(exp.norm, na.rm=T), .N),
                           list(var=var, cds_len=round(log10(cds_len), 1))]


### plot
  genome_wide <- ggplot(data=genebody.ag, aes(x=var, y=cds_len, fill=mu)) +
  geom_tile() +
  scale_fill_viridis() +
  xlab("Distance along genebody 5'->3'") +
  ylab("log10(CDS Length)") +
  theme_bw()


  gene <- ggplot() +
  geom_ribbon(data=genebody.ag[cds_len==round(log10(12540), 1)], aes(x=var, ymax=mu-2*sd, ymin=mu+2*sd), alpha=.3) +
  geom_line(data=genebody.ag[cds_len==round(log10(12540), 1)], aes(x=var, y=mu)) +
  geom_line(data=genebody[GENE_ID=="Daphnia00787"], aes(x=var, y=exp.norm, group=samp, color=superclone)) +
  ylab("Normalized gene expression") +
  xlab("Distance along genebody 5'->3'") +
  theme_bw()

  genome_wide + gene 
