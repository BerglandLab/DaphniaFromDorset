module load  gcc/7.1.0  openmpi/3.1.4
module load R/3.5.3
R

#!/usr/bin/env Rscript

    ### libraries
      library(data.table)
      library(foreach)
      library(ggplot2)
      library(tidyverse)
      library(doMC)
      registerDoMC(20)

### load files

  load("totalgenomewidenoognoACsub.Rdata")

### Graph
  IBSplot <- ggplot(data=totalgenomewidenoognoACsub, aes(x=type, y=mean_ratiosim, group=type)) + geom_boxplot() +
    geom_hline(yintercept=1, color="red") + theme_bw() +
    ylab("AvsC IBS/Comparison IBS") + xlab("Comparison")

  IBSplotB <- IBSplot + scale_x_discrete(labels=c("Within D8", "Within DCat/DBunk", "Between DCat/D8/DBunk", "Kilwood vs D10", "Kilwood vs W1", "Kilwood vs W6"))

  IBSgenomewide <- IBSplotB + coord_flip()

  ggsave(IBSgenomewide, file="IBSgenomewide.pdf")
