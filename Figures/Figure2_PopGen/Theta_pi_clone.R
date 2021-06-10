# Connor Murray
# 7.15.2020

# Libraries
library(data.table)
library(tidyverse)
library(foreach)
library(cowplot)

#setwd("C:/Users/Conno/Desktop/Spring2020/rohan/data/")
setwd("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure2")

# Filter metadata sample file
samps <- fread("cloneinfo.txt")
samps <- data.table(samps %>% filter(Nonindependent==0 & is.na(LowReadDepth) & Species=="pulex"))

# Summary output from ROHan
out <- readRDS("summary.all.rds")
out[pond=="Dcat"]$pond <-"DCat"

# Subsample 1 clone per SC
setkey(1234)
rand.sc <- data.table(samps %>% filter(population %in% c("D8","DBunk","DCat")) %>%
                        group_by(SC) %>% sample_n(1))$clone

# OO clones of interest
oo.sc <- data.table(samps %>% filter(population %in% c("D8","DBunk","DCat") & SC %in% "OO"))$clone

out$pond <- factor(out$pond, levels = c("DBunk", "D8", "DCat"))


# Parent-offspring relationships
par.off <- data.table(rbind(read.csv("D8ParentOffspringRelationships.long.csv"),
                            read.csv("DBunkParentOffspringRelationships.long.csv")))

# A x C cross clones
cross <- data.table(read.csv("A.C.cross.csv", header = TRUE))


# Selfing parents
selfing.par <- c(samps[SC %in% c("poB", "poH", "poC", "poW")]$clone)

# Add column for graph
out[clone %in% selfing.par, SC.share:=str_replace(paste(SC), "po", "")]
out[!clone %in% selfing.par, SC.share:=paste(SC)]

### save output
save(out, rand.sc, oo.sc, cross, selfing.par, file="genome_theta.Rdata")






























# Cowplot of all graphs
plot_grid(a, plot_grid(b, c, labels = c("B.", "C."),
                       rel_widths = c(1,2), align="h", axis = "tb"),
          nrow = 2, labels = c("A.", ""))















###roh plots

a <- ggplot(out[clone %in% c(rand.sc, oo.sc)][data %in% c("theta_inc_ROH", "seg_in_ROH%")],
       aes(x=year, y=global)) +
  geom_boxplot() +
  facet_grid(data~pond, scales="free_y") +
  labs(x="Year", y="Theta pi") +
  theme_classic() +
  theme(axis.text.x =element_text(face="bold", size=12),
        axis.text.y = element_text(face="bold", size=12),
        axis.title.x = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=15),
        title = element_text(face="bold", size=15),
        strip.text = element_text(face="bold", size=12))

        b <- ggplot() +
          geom_point(data=out[!clone %in% (cross$clone)][SC %in% c("A", "C")][data %in% c("theta_inc_ROH", "seg_in_ROH%")],
                     aes(x=SC, y=global, color=SC), size=3) +
          geom_point(data=out[clone %in% cross[cross=="F1"]$clone][!SC=="AL"][data %in% c("theta_inc_ROH", "seg_in_ROH%")],
                     aes(x="F1-AxC", y=global), color="purple", size=4) +
          geom_point(data=out[clone %in% cross[cross=="F2"]$clone][data %in% c("theta_inc_ROH", "seg_in_ROH%")],
                     aes(x="F2-AxC", y=global), color="orange", size=4) +
          labs(x="", y="Theta pi") +
          theme_bw() +
          facet_grid(data~., scales="free_y") +

          theme(legend.position = "none",
                axis.text.x =element_text(face="bold", size=12),
                axis.text.y = element_text(face="bold", size=12),
                axis.title.x = element_text(face="bold", size=15),
                axis.title.y = element_text(face="bold", size=15))
