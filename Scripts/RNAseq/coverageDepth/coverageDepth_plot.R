#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

library(data.table)
library(foreach)

fn <- list.files("/scratch/aob2x/daphnia_hwe_sims/rnaseq/coverage/", "coveragePos.delim", full.names=T)

dat <- foreach(fn.i=fn, .combine="rbind")%do%fread(fn.i)
save(dat, file="~/coverage_pos_rna.Rdata")


scp aob2x@rivanna.hpc.virginia.edu:~/coverage_pos_rna.Rdata ~/.

library(data.table)
library(ggplot2)

load("~/coverage_pos_rna.Rdata")
setnames(dat, names(dat), c("samp", "chrLen", "readsMapped", "chr", "pos", "rd"))
dat[,rd.norm:=rd/readsMapped]
dat[pos>=5191562 & pos<=5204101, gene:="Daphnia00787"]


ggplot(data=dat[], aes(x=pos, y=rd.norm, group=samp, color=samp)) +
geom_line() +
geom_hline(yintercept=0)


plot(V4~V3, dat)
abline(h=0)
