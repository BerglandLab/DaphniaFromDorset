#!/usr/bin/env Rscript

### libraries
  library(data.table)
  library(ggplot2)

### Load file
  shannonh <- fread("ShannonsHforR.csv")

### Graph
  ShannonH <- ggplot(data=shannonh, aes(x=Year, y=ShannonsH, group=Pond, color=Pond)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

    ggsave(ShannonH, file="ShannonH.pdf")
