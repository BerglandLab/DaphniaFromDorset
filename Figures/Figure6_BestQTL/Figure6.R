### libraries
  library(ggplot2)
  library(data.table)
  library(cowplot)
  library(patchwork)
  library(ggthemes)
  library(foreach)
  library(DESeq2)

### load data
  setwd("/Users/alanbergland/Documents/GitHub/")
  load(file="DaphniaPulex20162017Sequencing/AlanFigures/Figure4/qtl_polar_withNames.Rdata")
  load("DaphniaPulex20162017Sequencing/AlanFigures/Figure4/f1_pool_polar.Rdata")
  load(file="DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/differential_expression_supplement.Rdata")
  load(file="DaphniaPulex20162017Sequencing/AlanFigures/Figure4/differentialExpression_withNames.Rdata")

  qtlcoding <- qtl.polar[,list(unique(PA42qtl), qtl=unique(qtl)), list(final_QTL_ID)]



##############################
### some data manipulation ###
##############################

  ### pond level difference
      ab <- qtl.polar[,list(geno=c("male_male", "male_pe", "pe_pe")[which.max(c(n.male_male, n.male_pe, n.pe_pe))],
                            size=N[1]), list(clone, pond=pond.y, qtl, PA42qtl, sc=SC.uniq, year=year, final_QTL_ID)]

      abr <- ab[,list(geno=rep(geno, size)), list(clone, pond, sc, qtl, PA42qtl, year, final_QTL_ID)]

      abrf <- abr[,list(male_freq=(2*sum(geno=="male_male") + sum(geno=="male_pe"))/(2*length(geno)), n=2*length(geno)),
                    list(pond, year, qtl, PA42qtl, final_QTL_ID)]
      abrf[,se:=male_freq*(1-male_freq)/sqrt(n)]

      qtl12.pond <- abrf[pond%in%c("D8", "DBUNK", "DCAT")][year>2016 & final_QTL_ID==12]
      qtl12.pond[pond=="DBUNK", pond:="DBunk"]
      qtl12.pond[pond=="DCAT", pond:="DCat"]

      qtl12.pond[,pond:=factor(pond, levels=c("DBunk", "D8", "DCat"))]

  ### F1 data
    f1.pool.merge <- merge(f1.pool.merge, qtlcoding, by="qtl")

    ### figure version
    m.ag <- f1.pool.merge[,list(propmale=sum(propmale*N)/sum(N), sd=sd(propmale), nObs=sum(N), nClone=length(clone)),
                           list(qtl, geno, final_QTL_ID)]

    m.ag[,se:=(propmale*(1-propmale))/sqrt(nObs)]
    m.ag[,lci:=propmale-1.96*se]
    m.ag[,uci:=propmale+1.96*se]

    m.ag[geno=="male_male", geno:="male+/male+"]
    m.ag[geno=="male_pe", geno:="male+/male-"]
    m.ag[geno=="pe_pe", geno:="male-/male-"]

    ### linear model
    f1.pool.merge[final_QTL_ID==12]
    t1 <- glmer(propmale~geno+gr+(1|clone), data=f1.pool.merge[final_QTL_ID==12],
              family=binomial(), weights=f1.pool.merge[final_QTL_ID==12]$N)
    t2 <- glmer(propmale~gr+(1|clone), data=f1.pool.merge[final_QTL_ID==12],
              family=binomial(), weights=f1.pool.merge[final_QTL_ID==12]$N)
    anova(t1, t2)
)
  ### gene counts

    gene_counts <- foreach(i=dec[qLFC.noCN.goodChr>.9 & !is.na(old_QTL_ID)]$GeneID, .combine="rbind")%do%{
      #i<-"Daphnia00787"
      gene_counts1 <- as.data.table(plotCounts(dds, gene=which(rownames(dds)==i), returnData=T, intgroup="clone", transform=T, normalized=T))
      gene_counts2 <- as.data.table(plotCounts(dds, gene=which(rownames(dds)==i), returnData=T, intgroup="superclone", transform=T, normalized=T))
      gene_counts <- merge(gene_counts1, gene_counts2, by="count")
      gene_counts[,gene:=i]
    }
    gene_counts <- merge(gene_counts, dec, by.x="gene", by.y="GeneID")





############
### plot ###
############

  qtl12.pond.plot <-
  ggplot(data=qtl12.pond,
          aes(x=as.factor(year), y=male_freq, group=pond, color=pond)) +
  geom_errorbar(aes(ymin=male_freq-2*se, ymax=male_freq+2*se), width=0.05) +
  geom_point(size=2.5)+
  geom_line() +
  xlab("Year") + ylab("Freq male allele") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
  scale_color_colorblind()


  qtl12.f1s <-
  ggplot(data=m.ag[final_QTL_ID==12],
          aes(x=geno, y=propmale, group=qtl, linetype=as.factor(qtl))) +
  geom_line() +
  geom_point(size=2.5) +
  geom_errorbar(aes(ymin=propmale-2*se, ymax=propmale+2*se), width=0.2) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Genotype") + ylab("Proportion male") +
  geom_linerange(aes(ymin=lci, ymax=uci))

  qtl_12_de <-
    ggplot(gene_counts[final_QTL_ID==12][gene=="Daphnia00787"],
            aes(x=superclone, y=count, color=superclone, shape=clone)) +
    geom_point(size=3) +
    ylab("Normalized expression") +
    facet_grid(.~gene, scales="free") + theme_bw()



#### mega plot
layout <- "
AB
CD
"


qtl.confirm <-
qtl12.pond.plot + qtl12.f1s +
qtl_12_de + plot_spacer() +
plot_layout(design = layout) +
plot_annotation(tag_levels = 'A')

qtl.confirm


ggsave(qtl.confirm, file="~/qtl_confirm.pdf")

























library(GenomicAlignments)
library(ggbio)
ga <- readGAlignments("/Users/alanbergland/bam_slices/D86A.forward.bam",
                      use.names=TRUE,
                      param=ScanBamParam(which=GRanges("Scaffold_2217_HRSCAF_2652", IRanges(5188282, 5215770), flag=flag0)))
cov <- coverage(unlist(ga))


cvg2 <- coverage(readGAlignments(f1, param=ScanBamParam(flag=flag0))) stopifnot(identical(cvg, cvg2))


extractCoverageFromBAM <- function(bamfile) {
  stopifnot(is(bamfile, "BamFile"))
  ## This ScanBamParam object allows us to load only the necessary
  ## information from the file.
  param <- ScanBamParam(flag=flag0, what=c("rname", "pos", "cigar"))
  bam <- scanBam(bamfile, param=param)[[1]]
  ## Note that unmapped reads and reads that are PCR/optical duplicates ## have already been filtered out by using the ScanBamParam object
  ## above.
  f <- factor(bam$rname, levels=seqlevels(bamfile))
  irl <- extractAlignmentRangesOnReference(bam$cigar, pos=bam$pos, f=f)
  coverage(irl, width=seqlengths(bamfile))
}

extractCoverageFromBAM("/Users/alanbergland/bam_slices/D86A.forward.bam")

greads <- readGappedReads("/Users/alanbergland/bam_slices/D86A.forward.bam")
greads
tmp <- qseq(greads)






library(Rsamtools)
library(GenomicAlignments)
bfl <- BamFile(system.file(package="Rsamtools", "extdata", "ex1.bam"))
param1 <- ScanBamParam(which=GRanges("seq1", IRanges(10, width=1)))
param2 <- ScanBamParam(which=GRanges("seq1", IRanges(c(10, 10), width=1)))




bfl <- BamFile("/Users/alanbergland/bam_slices/D86A.reverse.bam")
pu <- pileup(bfl,
            ScanBamParam(which=GRanges("Scaffold_2217_HRSCAF_2652", IRanges(5188282, 5215770))))
pu <- as.data.table(pu)
ggplot(data=pu, aes(x=pos, y=count, group=strand, color=strand)) + geom_line()




### gene expression plot
  txdb <- makeTxDbFromGFF(file="/Users/alanbergland/daphnia_ref/Daphnia.aed.0.6.gff", format="gff3")

  ### forward strand
    ga <- readGAlignments("/Users/alanbergland/bam_slices/D86A.forward.bam",
                          use.names=TRUE,
                          param=ScanBamParam(which=GRanges("Scaffold_2217_HRSCAF_2652", IRanges(5188282, 5215770))))

    p1 <- autoplot(ga, geom = "rect")
    p2 <- autoplot(ga, geom = "line", stat = "coverage")

    p3<- autoplot(txdb,
               which=GRanges("Scaffold_2217_HRSCAF_2652", IRanges(5188282, 5215770), strand="+"),
               names.expr = "gene_id")



    forward.plot <- tracks(Reads=p1, Coverage=p2, Transcripts=p4, heights = c(0.3, 0.2, 0.1)) + ylab("")

  ### reverse strand
    ga <- readGAlignments("/Users/alanbergland/bam_slices/D86A.reverse.bam",
                          use.names=TRUE,
                          param=ScanBamParam(which=GRanges("Scaffold_2217_HRSCAF_2652", IRanges(5188282, 5215770))))

    p1 <- autoplot(ga, geom = "rect")
    p2 <- autoplot(ga, geom = "gapped.pair", stat = "coverage")

    p3<- autoplot(txdb,
               which=GRanges("Scaffold_2217_HRSCAF_2652", IRanges(5188282, 5215770), strand="-"),
               names.expr = "gene_id")



    reverse.plot <- tracks(Reads=p1, Coverage=p2, Transcripts=p4, heights = c(0.3, 0.2, 0.1)) + ylab("")



forward.plot / reverse.plot

pos <- sapply(coverage(gr[strand(gr)=="+"]), as.numeric)
pos <- data.frame(Chr=rep(names(pos), sapply(pos, length)), Strand=rep("+", length(unlist(pos))), Position=unlist(sapply(pos, function(x) 1:length(x))), Coverage=as.numeric(unlist(pos)))
neg <- sapply(coverage(gr[strand(gr)=="-"]), as.numeric)
neg <- data.frame(Chr=rep(names(neg), sapply(neg, length)), Strand=rep("-", length(unlist(neg))), Position=unlist(sapply(neg, function(x) 1:length(x))), Coverage=-as.numeric(unlist(neg)))
covdf <- rbind(pos, neg)
p <- ggplot(covdf, aes(Position, Coverage, fill=Strand)) +
	    geom_bar(stat="identity", position="identity") + facet_wrap(~Chr)
p



p2 <- autoplot(ga, geom = "line", stat = "coverage")
vcf <- readVcf(file="data/varianttools_gnsap.vcf", genome="ATH1")
p3 <- autoplot(vcf[seqnames(vcf)=="Chr5"], type = "fixed") + xlim(4000, 8000) + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y=element_blank())
txdb <- makeTxDbFromGFF(file="./data/TAIR10_GFF3_trunc.gff", format="gff3")
p4 <- autoplot(txdb, which=GRanges("Chr5", IRanges(4000, 8000)), names.expr = "gene_id")
tracks(Reads=p1, Coverage=p2, Variant=p3, Transcripts=p4, heights = c(0.3, 0.2, 0.1, 0.35)) + ylab("")
