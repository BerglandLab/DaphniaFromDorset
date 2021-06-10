#############################
#For Alan -- how to make the hapnet of 787
#############################

#Load libraries to R
library(pegas)
library(adegenet)
library(tidyverse)
library(magrittr)
library(factoextra)
library(FactoMineR)

#############################
# Load the R object
#############################
setwd("/Users/alanbergland/Documents/GitHub/DaphniaPulex20162017Sequencing/AlanFigures/Figure5")
load("dataForHapNet.Rdata")

#############################
# Remove obtusa and pulicaria
#############################
#DNA787 = DNA787[-grep("pulicaria",rownames(DNA787)),]
#DNA787 = DNA787[-grep("obtusa",rownames(DNA787)),]

#############################
#re format as data frame
#############################
tab(DNA787df) %>%
  as.data.frame() -> DNA787df

#############################
#add sample names to df
#############################
DNA787df %<>%
  mutate(clone =rownames(.),
         full_name = rownames(.))

#############################
# clean up names for prettiness
#############################
DNA787df$clone = gsub("hap[01]_","",DNA787df$clone)
DNA787df$clone = gsub("_787","",DNA787df$clone)

#############################
#Merge DF with metadata
#############################
left_join(DNA787df, metadata) %>%
  .[complete.cases(.$population),] -> DNA787df_mdt

#############################
#Construct the haplonet
#############################
table(rownames(DNA787))
DNA787Haps <- haplotype(DNA787)
DNA787Haps ## <--- this will print the haps to screen

#############################
#Run the network algorithm here:
#############################
DNA787Net <- haploNet(DNA787Haps)

#############################
#Concert the metadata into color schemes for the pie graphs!
#This is done multiple times for various metrics....
#############################
utils::stack(setNames(attr(DNA787Haps, "index"), rownames(DNA787Haps))) -> ind_haps

ind_haps %<>% mutate(hap_name = paste("hap",ind, sep = "_"))
ind_haps %<>% mutate(ind_name = row.names(DNA787)[.$values])

ind_haps ## <--- this is your base file with the haps of each individual

#############################
#pie colors for netwrork at the individual level:
#############################
ind.hap<-with(
  utils::stack(setNames(attr(DNA787Haps, "index"), rownames(DNA787Haps))),
  table(hap=ind, pop=rownames(DNA787)[values])
)

#############################
## This generate colors at the population level
#############################
pop.peakState<-with(
  utils::stack(setNames(attr(DNA787Haps, "index"), rownames(DNA787Haps))),
  table(hap=ind, pop=as.factor(DNA787df_mdt$"5881.A/G"[ind_haps$values]))
)

#############################
## This generate colors at the pond level
#############################
pop.pond<-with(
  utils::stack(setNames(attr(DNA787Haps, "index"), rownames(DNA787Haps))),
  table(hap=ind, pop=as.factor(DNA787df_mdt$population[ind_haps$values]))
)

#############################
# Now lets do some plotting #
#############################

#############################
#Plot by pond
#############################

pdf("~/DNA787Net.ponds.pdf") ## <---- save file begin
plot(DNA787Net,
     size = attr(DNA787Net, "freq"),
     fast = F,
     show.mutation = 0,
     cex = 0.1,
     pie = pop.pond,
     labels= F,
     scale.ratio = 0.6,
     threshold = 0)
legend("bottomright",
       sort(unique(DNA787df_mdt$population)),
       fill = rainbow(length(unique(DNA787df_mdt$population))),
       cex = 0.9)
dev.off() ## <---- save file end


,
col = rainbow(length(unique(DNA787df_mdt$population)))

#############################
#Plot by pop.peakState, i.e. mutation state
#############################
pdf("~/DNA787Net.pop.peakState.pdf") ## <---- save file begin
plot(DNA787Net,
size = attr(DNA787Net, "freq"),
fast = F,
show.mutation =0,
cex = 0.4,
pie=pop.peakState,
labels= F,
scale.ratio = 0.6,
threshold = 0)




,
col = rainbow(length(unique(DNA787df_mdt$"5881.A/G"))))
legend("bottomright",
       as.character(sort(unique(DNA787df_mdt$"5881.A/G"))),
       fill = rainbow(length(unique(DNA787df_mdt$"5881.A/G"))),
       cex = 0.9)
dev.off()  ## <---- save file end

#############################
#Plot by pop.peakState + mutations steps!!
#This one is a messy plot
#############################
pdf("DNA787Net.pop.peakState.muts.pdf") ## <---- save file begin
plot(DNA787Net,
     size = attr(DNA787Net, "freq"),
     fast = FALSE,
     show.mutation =3,
     cex = 0.4,
     pie=pop.peakState,
     labels= F,
     scale.ratio = 0.6,
     threshold = 0,
     col = rainbow(length(unique(DNA787df_mdt$"5881.A/G"))))
dev.off() ## <---- save file end




pdf("~/haplotype.pdf")
plot(DNA787Net,
size = attr(DNA787Net, "freq"),
fast = F,
show.mutation =0,
cex = 1,
labels= F,
scale.ratio = 0.6,
threshold = 0)
dev.off()
