library(data.table)
library(foreach)
library(ggplot2)

dat <- fread("~/allGenes.results.txt.gz")
dat[geneID=="Daphnia00787"]
