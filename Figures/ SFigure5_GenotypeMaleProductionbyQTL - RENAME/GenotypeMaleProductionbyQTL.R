### libraries
  library(ggplot2)
  library(data.table)
  library(foreach)
  library(patchwork)
  library(lme4)

setwd("/Users/alanbergland/Documents/GitHub/")

### Load in PA42 chr key
    PA42 <- fread("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/5kbchrassignHiCnew.csv")

### Load in qtl_polar (wild daps only: `DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_superclone_test/analysis.polarizeHapolotypes_poolSeq.R`)
    load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/f1_pool_polar.Rdata")

### load F1 cross data (made by `DaphniaPulex20162017Sequencing/AlanFigures/Figure4/makeData.F1_pheno.R`)
    load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/F1_pheno.Rdata")


############################
### male proportion      ###
############################

  male$gr <- ifelse(male$SC=="selfedA", "A", ifelse(male$SC=="B", "B", male$gr))
  male.ag <- male[!is.na(gr),list(propmale=sum(Males)/sum(NewTotal),
                        N=sum(NewTotal)),
                        list(clone, gr)]
  male.ag[,se:=sqrt((propmale*(1-propmale))/N)]

  male.ag[,lci:=propmale-1.96*se]
  male.ag[,uci:=propmale+1.96*se]

  male.ag$gr <- factor(male.ag$gr, levels=c("A", "AxC", "C", "CxC"))

#############################
### QTL polarization      ###
#############################
qtlcoding <- data.table(qtl=c(1:14), PA42qtl=c(3, 4, 2, 1, 9, 10, 11, 13, 14, 12, 8, 6, 7, 5))

f1.pool.mergeB <- merge(f1.pool.merge, qtlcoding, by="qtl")
f1.pool.merge <- f1.pool.mergeB

m.ag <- f1.pool.merge[,list(propmale=mean(propmale), sd=sd(propmale)), list(qtl, geno, PA42qtl, gr)]
sampsize <- length(unique(f1.pool.merge$clone))

m.ag[,se:=sd/(sqrt(sampsize))]
m.ag[,lci:=propmale-1.96*se]
m.ag[,uci:=propmale+1.96*se]

#m2.ag <- m.ag[,list(propmale=mean(propmale), sd=sd(propmale)), list(geno)]

f1.pool.mergesub <- f1.pool.merge[, c("clone", "PA42qtl", "geno"), with=TRUE]
malesqtl <- merge(f1.pool.mergesub, male, by="clone", allow.cartesian=TRUE)

anovas <- foreach(q=1:14, .combine="rbind")%do%{
    f1sub <- malesqtl[PA42qtl==q]
    f1subout_a <- glmer(propmale~gr + (1|clone) + (1|Replicate), data=f1sub, family=binomial(), weights=NewTotal)
    f1subout_b <- glmer(propmale~gr + geno+(1|clone) + (1|Replicate), data=f1sub, family=binomial(), weights=NewTotal)
    aovcompare <- anova(f1subout_a, f1subout_b)
    p <- aovcompare[[8]][2]
    tmp <- data.table(PA42qtl=q, p=p)
    tmp
    }

  m.agsig <- merge(m.ag, anovas, by="PA42qtl")

  m.agsig$genoB <- ifelse(m.agsig$geno=="pe_pe", "male-/male-", ifelse(m.agsig$geno=="male_pe",
    "male+/male-", "male+/male+"))


  F1polarnew <- ggplot(data=m.agsig[p <= 0.05], aes(x=genoB, y=propmale, group=gr, color=gr)) + geom_line() +
    theme_bw() + theme(legend.position = "none") + xlab("Genotype") + ylab("Proportion male") +
    geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
    facet_wrap(~PA42qtl)

  ggsave(F1polarnew, file="./DaphniaPulex20162017Sequencing/AlanFigures/SFigure12_GenotypeMaleProductionbyQTL/F1polarnew.pdf")
