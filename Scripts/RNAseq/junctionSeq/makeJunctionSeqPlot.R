#module load gcc; /home/aob2x/R-3.3.3/bin/R

library(JunctionSeq)

load("~/jscs.Rdata")


buildAllPlots(jscs=jscs,
  outfile.prefix = "~/Daphnia00787",
  use.plotting.device = "png", gene.list="Daphnia00787", writeHTMLresults=T
)
