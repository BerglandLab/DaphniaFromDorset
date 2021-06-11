# pN/pS plot
# 6.11.2021
# by Connor Murray

# Libraries
library(data.table)
library(tidyverse)

setwd("C:/Users/Conno/Desktop/pnps-figure")

# pnps data for non-OO clones
pnps.out <- data.table(read.csv("pnps.nonOO.csv"))

# Metadata
samps <- data.table(read.csv("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.csv"))
samps <- samps[PnPsAnalysis%in%c(1,2)]
pnps.out <- data.table(pnps.out %>% full_join(samps %>% select(SC, PnPsAnalysis, year) %>% 
                              group_by(SC) %>% 
                              summarize(PnPsAnalysis=unique(PnPsAnalysis), 
                                        year=unique(year))))

# Merge OO data
pnps.out <- data.table(pnps.out %>% full_join(data.table(read.csv("pnps.OO.csv"))))

# Change name of SC OO
pnps.out[SC=="OO"]$SC <- "Outbred"

# Change names of year 
pnps.out[year %in% c(2016, 2017, 17), year.an:="2016/2017"]
pnps.out[year %in% c(2018, 2019, 18), year.an:="2018/2019"]

# Number sampled per bootstrap
samp <- unique(pnps.out$sample.n)
mac <- samp*2

# Add pnps and filters
pnps.out[af > 0.5, maf:= 1-af]
pnps.out[is.na(maf), maf:=af]
pnps.out[, nMissing:= (samp*2) - MAC/af]
pnps.out <- pnps.out[nMissing==0]
pnps.out.b <- pnps.out[MAC > samp, maf:=maf]

# Recalculate pnps for folded SFS
pnps.out <- pnps.out.b[, list(pn=sum(pn), ps=sum(ps)), 
                     list(maf=as.numeric(as.character(maf)), SC, rep, PnPsAnalysis, year.an)]
pnps.out <- pnps.out[, pnps:=pn/ps]

# Removes infinte values and zeros
pnps.out <- pnps.out[!pnps > 2.5 & !pnps==0] 

# Aggregate pnps by SC
pnps.out.b <- pnps.out.b[, list(pn=sum(pn), ps=sum(ps)), 
                       list(af=as.numeric(as.character(af)), SC, rep)]
pnps.out.b <- pnps.out.b[, pnps:=pn/ps]
sfs.pnps.ag <- pnps.out.b[,list(pn=mean(pn, na.rm=T),
                              ps=mean(ps, na.rm=T),
                              pnps=mean(pnps, na.rm=T),
                              lci=quantile(pnps, 0.025, na.rm=T), 
                              uci=quantile(pnps, 0.975, na.rm=T)),
                         list(SC, af=af)]

# Wide to long format
sfs.pnps.ag <- data.table(sfs.pnps.ag %>% gather(measure, measurement, pn:uci, factor_key=TRUE))

# Round af
sfs.pnps.ag <- data.table(sfs.pnps.ag %>% mutate(af=round(af, 2)))

# Aggregate pnps by SC
fsfs.pnps.ag <- pnps.out[,list(pn = mean(pn, na.rm=T),
                          ps = mean(ps, na.rm = T),
                          pnps = mean(pnps, na.rm = T),
                          lci= quantile(pnps, .025, na.rm=T), 
                          uci= quantile(pnps, .975, na.rm=T)),
                    list(SC, maf=as.factor(maf), PnPsAnalysis, year.an)]

# Make maf numeric
fsfs.pnps.ag <- data.table(fsfs.pnps.ag %>% mutate(maf=as.numeric(as.character(maf))))

# wide to long format
fsfs.pnps.ag <- data.table(fsfs.pnps.ag %>% gather(measure, measurement, pn:uci, factor_key=TRUE))

# Round maf
fsfs.pnps.ag <- data.table(fsfs.pnps.ag %>% mutate(maf=round(maf,2)))

# Breaks for maf
breaks.bin = c(2/mac, (samp-2)/mac, (samp-1)/mac, samp/mac)
fsfs.pnps.ag1 <- fsfs.pnps.ag[measure=="pnps"]

# Bin maf and take mean and CI of pnps
bin <- data.table(aggregate(fsfs.pnps.ag1$measurement, 
                  by=list(cut(as.numeric(as.character(fsfs.pnps.ag1$maf)), 
                  breaks = breaks.bin, labels= NULL), 
                  fsfs.pnps.ag1$SC, fsfs.pnps.ag1$year.an), mean),
      aggregate(fsfs.pnps.ag1$measurement, by=list(cut(as.numeric(as.character(fsfs.pnps.ag1$maf)), 
                breaks = breaks.bin, labels= NULL), 
                fsfs.pnps.ag1$SC, fsfs.pnps.ag1$year.an), quantile, probs=0.025)$x,
      aggregate(fsfs.pnps.ag1$measurement, by=list(cut(as.numeric(as.character(fsfs.pnps.ag1$maf)), 
                breaks = breaks.bin, labels= NULL), 
                fsfs.pnps.ag1$SC, fsfs.pnps.ag1$year.an), quantile, probs=0.975)$x) 

# Renaming columns and rbinding 1/2N pnps
colnames(bin) <- c("maf", "SC", "year.an", "pnps", "lci", "uci")

# Reordering for ggplot
bin <- fsfs.pnps.ag[!measure%in%c("pn","ps")][maf>0][maf==min(maf), -3] %>% 
                        pivot_wider(names_from=measure, values_from=measurement) %>% 
                        rbind(bin)

# Change MAF names for facets
bin[maf %in% c(round(1/mac, 2))]$maf <- "1/2N"
bin[maf %in% c(paste("(",round((samp-1)/mac,3),",",samp/mac,"]",sep=""))]$maf <- "~ 0.5"

# Values for statistical comparisons
pnps.out <- data.table(pnps.out %>% mutate(maf.round=round(maf, 4)))

# Color panel
# show_col(hue_pal()(7))  

# Final pnps plot
ggplot(bin[maf %in% c("1/2N","~ 0.5")],
            aes(x=reorder(SC, desc(SC)), y=pnps, 
                fill=maf, group=interaction(SC,maf))) +
  geom_col(position =  position_dodge(1)) +
  geom_errorbar(aes(ymin=bin[maf %in% c("1/2N","~ 0.5")]$lci, 
                    ymax=bin[maf %in% c("1/2N","~ 0.5")]$uci),
                width=0.3, position = position_dodge(1)) +
  facet_wrap(~year.an, ncol = 1, scales = "free_y") +
  theme_classic() +
  theme(panel.background = element_rect(colour = "black", size=1),
        axis.text.x =element_text(face="bold", size=14),
        axis.text.y = element_text(face="bold", size=14),
        axis.title.x = element_text(face="bold", size=16),
        axis.title.y = element_text(face="bold", size=16),
        strip.text = element_text(face="bold", size=16),
        legend.position = c(0.88, 0.2),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(face="bold", size=14),
        legend.title = element_text(face="bold", size=16),
        legend.background = element_blank()) +
  labs(x="Clonal lineage", y="pN/pS", 
       fill="MAF bin") + 
  scale_x_discrete(position = "bottom") +
  coord_flip()

