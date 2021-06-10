### Load libraries
library(data.table)
library(ggplot2)
library(tidyverse)
library(foreach)
library(vegan)

### Load superclone file

sc <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")

### Keep only wild caught, independent, pulex samples from D8, DBunk, DCat, and D10

scB <- sc[Nonindependent==0 & LabGenerated==0 & Species=="pulex"]
scC <- scB[population=="D8" | population=="DBunk" | population=="DCat" | population=="Dcat"]
scC$population <- str_replace(scC$population, "Dcat", "DCat")

### Get into correct format
scC$SCB <- ifelse(scC$SC=="OO", scC$clone, scC$SC)

### Aggregate by pond/year
scC.ag <- scC %>% count(population, year, SCB)
scC.ag$popyear <- paste(scC.ag$population, scC.ag$year, sep="_")
scC.agsub <- scC.ag[, c("popyear", "SCB", "n"), with=FALSE]
scC.agsubwide <- pivot_wider(scC.agsub, names_from=SCB, values_from=n)
scC.agsubwide[is.na(scC.agsubwide)] <- 0

scC.agsubN <- scC.agsub %>% group_by(popyear) %>% summarise(sumN=sum(n))


### Calculate Shannon's diversity
shannondiv <- diversity(scC.agsubwide[,2:158], index="shannon")
pondyear <- scC.agsubwide$popyear

shannondivdt <- data.table(popyear=scC.agsubwide$popyear, shannondiv=shannondiv)

simpsondiv <- diversity(scC.agsubwide[,2:158], index = "simpson")
pondyear <- scC.agsubwide$popyear

simpsondivdt <- data.table(popyear=scC.agsubwide$popyear, simpsondiv=simpsondiv)
