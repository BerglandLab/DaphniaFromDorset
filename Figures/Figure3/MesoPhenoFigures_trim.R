#!/usr/bin/env Rscript

### libraries
        library(data.table)
        library(ggplot2)
        library(tidyverse)
        library(patchwork)

### Load data files

  load("meso7AC.Rdata")
  load("dat.ag.Rdata")
  dat.ag$Trmt <- str_replace(dat.ag$Trmt,"C", "control")


### Figures
  week7ACtotalind <- ggplot(data=meso7AC, aes(x=Week, y=TotalPop, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("# Individuals") + theme(legend.position = "none")

  week7ACmomPEnormWeek <- ggplot(data=meso7AC, aes(x=Week, y=propPEtotalind, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Asexual Females") +
    theme(legend.position = "none")

  week7ACMomwEppnormWeek <- ggplot(data=meso7AC, aes(x=Week, y=propMomEpp, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Sexual Females") +
    theme(legend.position = "none")

  week7ACPropMalenormWeek <- ggplot(data=meso7AC, aes(x=Week, y=PropTotMale, color=SC, group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Prop Male") +
    theme(legend.position = "none")


  week7ACEmbWeekB <- ggplot(data=meso7AC, aes(x=Week, y=log10(EstimatedEmbryosTot), color=SC,
    group=CloneRep, linetype=Clone)) + geom_line() +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("log10(# Sexual Embryos)") +
    scale_y_continuous(breaks=c(0, 1, 2, 3), labels=c("0", "10", "100", "1000")) +
    theme(legend.position = "bottom") + guides(linetype = FALSE) + labs(color="Superclone")


  Methyl <- ggplot(data=dat.ag, aes(x=PrePost, group=gr, y=propmale, color=SC, linetype=Clone)) +
    geom_line(position=position_dodge(0.3)) +
    scale_size_manual(values = c(control=0.5,MF=1)) +
    geom_errorbar(aes(ymin=lci, ymax=uci), , linetype=1, width=.3, position=position_dodge(0.3)) +
    geom_point(data=dat.ag, aes(x=PrePost, group=gr, y=propmale, color=SC, shape=Trmt),
    position=position_dodge(0.3), size=2) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    xlab("Pre/Post Treatment") + ylab("Prop Male") +
    guides(color=FALSE) + theme(legend.position = c(0.21, 0.8)) + labs(shape="Treatment") +
    guides(linetype = FALSE)



### Total graph?


  total <- week7ACtotalind + week7ACmomPEnormWeek + week7ACMomwEppnormWeek + week7ACPropMalenormWeek +
    week7ACEmbWeekB + Methyl + plot_layout(nrow = 2, byrow = TRUE, guides = 'collect') +
    plot_annotation(tag_levels = 'A')

  upper <- week7ACtotalind + week7ACmomPEnormWeek + week7ACMomwEppnormWeek +
    plot_annotation(tag_levels = 'A')

  lower <- week7ACPropMalenormWeek + week7ACEmbWeekB + Methyl +
    plot_annotation(tag_levels = 'A')

  totalB <- upper/lower + plot_layout(nrow = 2, byrow = TRUE) +
    plot_annotation(tag_levels = 'A')

  ggsave(totalB, file="Figure3.pdf")
