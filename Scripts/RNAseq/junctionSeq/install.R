module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0

#Install CRAN packages:
install.packages("statmod")
install.packages("plotrix")
install.packages("stringr")
install.packages("Hmisc")
install.packages("locfit")
#Install Bioconductor packages:
source("http://bioconductor.org/biocLite.R");
biocLite();
library("Biobase");
library("BiocGenerics");
library("BiocParallel");
library("GenomicRanges");
library("IRanges");
library("S4Vectors");
library("genefilter");
library("geneplotter");
library("SummarizedExperiment");
library("DESeq2");


install.packages("http://hartleys.github.io/JunctionSeq/install/JunctionSeq_LATEST.tar.gz",
                   repos = NULL,
                   type="source")
