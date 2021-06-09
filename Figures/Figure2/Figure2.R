### libraries
  library(data.table)
  library(ggplot2)
  library(viridis)
  library(tidyverse)
  library(ggrepel)
  library(cowplot)
  library(patchwork)
  library(ggprism)

### set working directory
  setwd("/Users/alanbergland/Documents/GitHub")

############
### data ###
############

### load IBS data
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure2/ibs.long.Rdata") ### loads `ibs.long` object

### load diversity estimates (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure2/Theta_pi_clone.R`)
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure2/genome_theta.Rdata")

### load ROH information (made by )
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure2/rohan_summary.Rdata")


#######################
### plot components ###
#######################

### IBS
  h.just <- .25
  v.just <- .25
  l.size <- 1.5

  corrmatrix <- ggplot(data=ibs.long, aes(scid.a, scid.b, fill=IBS)) +
                geom_raster() +
                scale_fill_viridis(option="D") +
                theme(axis.title=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank())

### Theta-pi
  ## among ponds
    theta_ponds <-
      ggplot(out[clone %in% c(rand.sc, oo.sc)][data %in% c("theta_inc_ROH")],
           aes(x=year, y=global)) +
      geom_boxplot() +
      facet_wrap(~factor(pond, levels=c("DCat", "D8", "DBunk"))) +
      xlab("Year") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
            axis.text.y = element_text(size=12),
            axis.title = element_text(size=14)) +
      scale_y_continuous(breaks=c(.002, .0025, .003, .0035), labels=c(2, 2.5, 3, 3.5), limits=c(0.0019, 0.00375)) +
      ylab(expression(paste(theta[pi], " (x", 10^-3, ")", sep="")))


  # AxC crosses
    theta_AC <-
    ggplot() +
    geom_point(data=out[!clone %in% (cross$clone)][SC %in% c("A", "C")][data %in% c("theta_inc_ROH")][pond=="D8"],
               aes(x=SC, y=global), size=2) +
    geom_point(data=out[clone %in% cross[cross=="F1"]$clone][!SC=="AL"][data %in% c("theta_inc_ROH")][pond=="D8"],
               aes(x="F1-AxC", y=global), size=2) +
    geom_point(data=out[clone %in% cross[cross=="F2"]$clone][data %in% c("theta_inc_ROH")][pond=="D8"],
               aes(x="F2-AxC", y=global), size=2) +
    theme_bw() +
    facet_wrap(~pond) +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size=12),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=14),
          legend.position = "none") +
    scale_y_continuous(breaks=c(.002, .0025, .003, .0035), labels=c(2, 2.5, 3, 3.5), limits=c(0.0019, 0.00375)) +
    ylab(expression(paste(theta[pi], " (x", 10^-3, ")", sep=""))) +
    xlab("")

# Selfed vs offspring (Supp Fig)
  theta_selfed <-
    ggplot() +
    geom_point(data=out[clone %in% selfing.par][data %in% c("theta_inc_ROH")],
               aes(x=as.factor(c("1", "1", "1", "1")),
                   y=global), size=2) +
    geom_point(data=out[SC %in% c("B", "H", "C", "W")][data %in% c("theta_inc_ROH")],
               aes(x=as.factor(SC), y=global), size=2) +
    labs(x="", y="") +
    scale_x_discrete(labels=c("B"="Offspring", "C"="Offspring",
                              "H"="Offspring", "W"="Offspring",
                              "1"="Parent")) +
    facet_grid(~SC.share, scales = "free_x") +
    guides(color=FALSE) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size=12, angle = 45, hjust=1, ),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=14)) +
    scale_y_continuous(breaks=c(.002, .0025, .003, .0035), labels=c(2, 2.5, 3, 3.5), limits=c(0.0019, 0.00375)) +
    ylab(expression(paste(theta[pi], " (x", 10^-3, ")", sep="")))



### ROHs
  rl.ag <- rl[ROH_LENGTH>=100000, list(nroh=.N, sroh=sum(ROH_LENGTH), mroh=mean(ROH_LENGTH)), list(clone, SC.uniq, population)]
  rl.ag.ag <- rl.ag[, list(nroh=mean(nroh), sroh=mean(sroh), mroh=mean(mroh)), list(SC.uniq, population)]


  rl.ag.ag[SC.uniq%in%c("A", "C") & population=="D8", lab:=SC.uniq]

  rl.ag.ag[is.na(lab), lab:=""]
  rl.ag.ag[population=="Dcat", population:="DCat"]
  #rl.ag.ag <- merge(rl.ag.ag, m.ag, by="SC.uniq")

  m.ag <- m[,list(.N), list(year, SC.uniq)]
  rl.ag.ag <- merge(rl.ag.ag, m.ag, by="SC.uniq")
  rl.ag.ag[year!="2017", lab:=""]

  ### between populations
    summary(lm(sroh~population, rl.ag.ag[population%in%c("D8", "DBunk")]))
    summary(lm(nroh~population, rl.ag.ag[population%in%c("D8", "DBunk")]))

  ### between years
    summary(lm(sroh~as.factor(year), rl.ag.ag[population%in%c("D8")]))
    summary(lm(nroh~as.factor(year), rl.ag.ag[population%in%c("D8")]))


    t1 <- lm(nroh~sroh, rl.ag.ag)

    rl.summary <- rl.ag.ag[,list(sroh.ave=mean(sroh), sroh.se=sd(sroh)/sqrt(.N),
                               nroh.ave=mean(nroh), nroh.se=sd(nroh)/sqrt(.N)),
                          list(year, population, lab)][population%in%c("D8", "DBunk")]




  roh_plot <-
  ggplot(data=rl.ag.ag[population%in%c("D8", "DBunk", "DCat")],
        aes(x=sroh, y=nroh, color=as.factor(year), label=lab, group=as.factor(year))) +
  geom_abline(aes(slope=coef(t1)[2], intercept=coef(t1)[1]))  +
  geom_point(size=1.5) +
  geom_segment(data=rl.summary[population%in%c("D8", "DBunk", "DCat")],
            aes(x=sroh.ave-1.96*sroh.se, xend=sroh.ave+1.96*sroh.se,
                y=nroh.ave, yend=nroh.ave,
                group=as.factor(year)),
            size=1, color="black") +
  geom_segment(data=rl.summary[population%in%c("D8", "DBunk", "DCat")],
            aes(x=sroh.ave, xend=sroh.ave,
                y=nroh.ave-1.96*nroh.se, yend=nroh.ave+1.96*nroh.se,
                group=as.factor(year)),
            size=1, color="black") +
  geom_point(data=rl.summary[population%in%c("D8", "DBunk", "DCat")],
            aes(x=sroh.ave, y=nroh.ave, fill=as.factor(year), group=as.factor(year)),
            size=5, colour="black",pch=21) +
  geom_label_repel(box.padding=1.5, color="black", size=5) +
  facet_wrap(~factor(population, levels=c("DCat", "D8", "DBunk"))) +
  theme_bw() +
  ylab("Number ROH") +
  xlab(expression(paste("Summed length ROH (x", 10^6, ")", sep=""))) +
  theme(legend.position = "none",
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title = element_text(size=24),
        strip.text = element_text(size = 24)) +
  scale_x_continuous(breaks=c(0, 2.5e6, 5e6, 7.5e6), labels=c(0, 2.5, 5, 7.5))





### pedigree
  pedigree_plot <-
  ggplot() +
  draw_image("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure2/Pedigree.png") +
  theme(axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())



#################
### MEGA-plot ###
#################


### trual
  layout <- "
  AAAAABBBB
  AAAAACCDD
  EEEEEEEEE
  FFFFFFFFF
  FFFFFFFFF"

  mega <-
  corrmatrix  + theta_ponds + theta_AC + theta_selfed+ roh_plot +  pedigree_plot +
  theme(strip.text=element_text(size=12), plot.tag = element_text(size=16)) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')

  ggsave(mega, file="~/mega_diversity_v3.pdf", height=11, w=8)

### this works

layout <- "
AABBEE
AACDEE
FFFFFF
FFFFFF
FFFFFF"

corrmatrix + theta_ponds + theta_AC + theta_selfed + roh_plot + pedigree_plot +
theme(strip.text=element_text(size=12), plot.tag = element_text(size=16)) +
plot_layout(design = layout) +
plot_annotation(tag_levels = 'A')


ggsave(mega, file="~/mega_diversity_v4.pdf")

ggsave(roh_plot + theme(strip.text = element_text(size = 18)), file="~/roh_plot.pdf", h=3, w=14)




###older
layout <- "
AABBC
AAEED
FFFFF
FFFFF
FFFFF"

layout <- "
AAFFFF
AAFFFF
BBBBEE
CCDDEE
"

mega <-
