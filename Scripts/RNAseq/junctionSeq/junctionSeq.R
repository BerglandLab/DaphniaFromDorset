# module load gcc/7.1.0  openmpi/3.1.4
#install.packages("BiocManager")
#BiocManager::install(version = "3.4")
#
#BiocManager::install("DESeq2", version = "3.4")

#BiocManager::install(version = "3.3")
#source("https://bioconductor.org/biocLite.R")
#biocLite("JunctionSeq")

#install.packages("https://cran.r-project.org/src/contrib/Archive/latticeExtra/latticeExtra_0.6-28.tar.gz", repos=NULL)
#install.packages("https://cran.r-project.org/src/contrib/Archive/Hmisc/Hmisc_4.0-3.tar.gz", repos=NULL)
#install.packages("https://cran.r-project.org/src/contrib/Archive/survival/survival_2.41-3.tar.gz", repos=NULL)
#install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.tar.gz", repos=NULL)
#install.packages("https://cran.r-project.org/src/contrib/Archive/plotrix/plotrix_3.6-2.tar.gz", repos=NULL)
#
#LoadR JunctionSeq:



  library("JunctionSeq")

#The sample decoder:
  decoder <- read.table("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/decoder.txt", header=T,stringsAsFactors=F)

#The count files:
  countFiles <- paste0("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/",
                      decoder$sample.ID,
                      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

#Run the analysis:
  jscs <- runJunctionSeqAnalyses(
    sample.files = countFiles,
    sample.names = decoder$sample.ID,
    condition = decoder$group.ID,
    use.novel.junctions=T,
    use.multigene.aggregates=T,
    flat.gff.file = "/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/withNovel.forJunctionSeq.gff.gz",
    verbose=TRUE, debug.mode = TRUE,
    nCores=20)

# save
  save(jscs, file="~/jscs.Rdata")

  # load(file="~/jscs.Rdata")
  # ijob -c5 -A berglandlab -p standard --mem 40
  writeCompleteResults(jscs, outfile.prefix="~/junctionSeq_out",
                      save.jscs = TRUE)
# plot
