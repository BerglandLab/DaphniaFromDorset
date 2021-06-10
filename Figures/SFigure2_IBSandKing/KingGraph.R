#!/usr/bin/env Rscript

### libraries
  library(data.table)
  library(foreach)
  library(ggplot2)
  library(tidyverse)
  library(patchwork)

### Load file
  load("kinshipscm2.Rdata")
  load("ibs.longunique.Rdata")

  kinshipscm2$type <- factor(kinshipscm2$type, levels=c("other", "clone", "selfedPO",
    "AxCF1vsACPO", "AxCF1fullsibs", "AvsC"))

### Graph
King <- ggplot(data=kinshipscm2[type=="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship)) + geom_point() +
  geom_point(data=kinshipscm2[type!="other" & medrdA > 14 & medrdB > 14], aes(x=IBS0, y=Kinship,
  color=factor(type, labels=c("Clonal", "Selfed parent-offspring",
  "Outcrossed parent-offspring", "Full-siblings", "A vs C")))) +
  labs(color="Type") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("#FF6666", "#CCCC00", "#00CC66", "#00CCFF", "#FF66FF")) +
  theme(legend.position = "none")

IBS <- ggplot(data=ibs.longunique, aes(x=IBS)) + geom_histogram(binwidth=0.001) +
  geom_vline(xintercept = 0.965, color="red") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

KingZoom <- ggplot(data=kinshipscm2[type=="other" & medrdA > 14 & medrdB > 14 &
  Kinship > 0 & IBS0 < 0.10], aes(x=IBS0, y=Kinship)) + geom_point() +
  geom_point(data=kinshipscm2[type!="other" & medrdA > 14 & medrdB > 14 &
  Kinship > 0 & IBS0 < 0.10], aes(x=IBS0, y=Kinship,
  color=factor(type, labels=c("Clonal", "Selfed parent-offspring",
  "Outcrossed parent-offspring", "Full-siblings")))) + scale_colour_discrete(drop = FALSE) +
  labs(color="Type") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position=c(0.7, 0.8)) +
  scale_color_manual(values=c("#FF6666", "#CCCC00", "#00CC66", "#00CCFF"))


kinshipscm2$clonallineage <- ifelse(kinshipscm2$SCcompare=="C_roC" |
  kinshipscm2$SCcompare=="C_C" | kinshipscm2$SCcompare=="roC_C" |
  kinshipscm2$SCcompare=="roC_roC", "C", ifelse(kinshipscm2$SCcompare=="B_roB" |
    kinshipscm2$SCcompare=="B_B" | kinshipscm2$SCcompare=="roB_B" |
    kinshipscm2$SCcompare=="roB_roB", "B", "other"))

exampleroC <- ggplot(data=kinshipscm2[clonallineage=="C" & medrdA > 10 & medrdB > 10],
  aes(x=IBS0, y=Kinship)) + geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylim(0, 0.5) + xlim(0,0.1) + theme(legend.position = "none")

exampleroB <- ggplot(data=kinshipscm2[clonallineage=="B" & medrdA > 10 & medrdB > 10],
  aes(x=IBS0, y=Kinship)) + geom_point() +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylim(0, 0.5) + xlim(0,0.1) + theme(legend.position = "none")


IBSplusKingB <- (IBS + KingZoom) / (exampleroC + exampleroB + King + plot_layout(widths=c(1,1,2))) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

ggsave(IBSplusKingB, file="IBSplusKingB.pdf")
