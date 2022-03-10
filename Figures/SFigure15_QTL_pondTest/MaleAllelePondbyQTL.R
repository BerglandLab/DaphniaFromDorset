### libraries
  library(ggplot2)
  library(data.table)
  library(foreach)
  library(patchwork)
  library(lme4)


### Load in PA42 chr key
    PA42 <- fread("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/5kbchrassignHiCnew.csv")

### Load in qtl_polar
    load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/qtl_polar.Rdata")

### Initial manipulation
  qtlcoding <- data.table(qtl=c(1:14), PA42qtl=c(3, 4, 2, 1, 9, 10, 11, 13, 14, 12, 8, 6, 7, 5))

  mqtl.polar <- merge(qtl.polar, qtlcoding, by="qtl")
  qtl.polar <- mqtl.polar
  setkey(mqtl.polar, PA42qtl)

  qtl.polar$geno <- ifelse(qtl.polar$n.male_male > qtl.polar$n.male_pe &
    qtl.polar$n.male_male > qtl.polar$n.pe_pe, "male+/male+",
    ifelse(qtl.polar$n.male_pe > qtl.polar$n.male_male &
    qtl.polar$n.male_pe > qtl.polar$n.pe_pe, "male+/male-",
    ifelse(qtl.polar$n.pe_pe > qtl.polar$n.male_male &
    qtl.polar$n.pe_pe > qtl.polar$n.male_pe,"male-/male-", "other")))

#### abundance genotype
      ab <- qtl.polar[,list(geno=c("male_male", "male_pe", "pe_pe")[which.max(c(n.male_male, n.male_pe, n.pe_pe))],
                            size=N[1]), list(clone, pond=pond.y, qtl, PA42qtl, sc=SC.uniq, year=year)]

      abr <- ab[,list(geno=rep(geno, size)), list(clone, pond, sc, qtl, PA42qtl, year)]

      abrf <- abr[,list(male_freq=(2*sum(geno=="male_male") + sum(geno=="male_pe"))/(2*length(geno)), n=2*length(geno)), list(pond, year, qtl, PA42qtl)]
      abrf[,se:=male_freq*(1-male_freq)/sqrt(n)]
      #abrf[,pond:=factor(pond, levels=c("DBUNK", "D8", "DCAT"))]

      abrf[,lci:=male_freq-1.96*se]
      abrf[,uci:=male_freq+1.96*se]

      good <- ggplot(data=abrf[pond%in%c("D8", "DBUNK", "DCAT")][year>2016], aes(x=as.factor(year), y=male_freq, group=pond, color=pond)) +
      geom_errorbar(aes(ymin=lci, ymax=uci), width=0.2) +
      geom_point() + geom_line() + facet_wrap(~PA42qtl) +
      xlab("Year") + ylab("Fequency of male+ allele") + theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

      ggsave(good, file="./DaphniaPulex20162017Sequencing/AlanFigures/SFigure12_GenotypeMaleProductionbyQTL/QTLbypond.pdf")
