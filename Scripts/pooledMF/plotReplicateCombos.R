# scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/gprime_peaks.replicates_combos.250K.05.Rdata ~/.

### libraries
  library(data.table)
  library(ggplot2)

  load("~/gprime_peaks.replicates_combos.250K.05.Rdata")

### plot
  ggplot(data=gprime, aes(x=POS, y=Gprime, color=CHROM)) +
  geom_line() +
  facet_grid(rep~CHROM)
