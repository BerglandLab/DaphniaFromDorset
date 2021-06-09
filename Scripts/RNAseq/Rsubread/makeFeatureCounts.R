# module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0; R

library(DESeq2)
library(Rsubread)
library(data.table)

### build feature counts table
  bamFiles <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam", "trim.bam", full.names=T)
  bamFiles <- bamFiles[!grepl("bam\\.", bamFiles)]
  saf <- flattenGTF("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf",
      method = "merge")

  fc <- featureCounts(files=bamFiles,
                annot.ext=saf,
                isPairedEnd=TRUE,
                nthreads=20,
                checkFragLength = FALSE)

  save(fc, file="/scratch/aob2x/daphnia_hwe_sims/rnaseq/featureCounts_trim.Rdata")
