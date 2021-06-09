## module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

#s
library(Rsubread)
#library(seqmagick)
#library(tidyverse)
#library(rtracklayer)
# Index pre-built during init.sh
# fa_read("/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa") %>% fa_write("/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.interleaved.fa", type = "interleaved")
# buildindex(basename="totalHiCwithallbestgapclosed.interleaved.rsubread.index", reference="/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.interleaved.fa", indexSplit=T)
#import("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff") %>% export("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf", "gtf")

args <- commandArgs(trailingOnly=TRUE)
filestem <- args[1]

#filestem="d8_179_1"

indexLoc <- "/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.interleaved.rsubread.index"

read1 <- paste("/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/", filestem, "_1.trim.fq.gz", sep="")
read2 <- paste("/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/", filestem, "_2.trim.fq.gz", sep="")
#read_merged <- paste("/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/", filestem, ".trim.merge.fq.gz", sep="")

saf <- flattenGTF("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf",
    method = "merge")

subjunc(index=indexLoc,
      readfile1=read1,
      readfile2=read2,
      readGroup = paste(filestem, ".pe", sep=""),
      input_format="gzFASTQ",
      output_format="BAM",
      output_file=paste("/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/", filestem, ".trim.bam", sep=""),
      nthreads=20,
      useAnnotation=T,
      annot.ext=saf)


#subjunc(index=indexLoc,
#      readfile1=read_merged,
#      readGroup = paste(filestem, ".merge", sep=""),
#      input_format="gzFASTQ",
#      output_format="BAM",
#      output_file=paste("/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/", filestem, "merge.bam", sep=""),
#      nthreads=20,
#      useAnnotation=T,
#      annot.ext=saf)
