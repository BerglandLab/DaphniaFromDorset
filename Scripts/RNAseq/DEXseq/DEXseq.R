


suppressPackageStartupMessages(library(DEXSeq))


#The sample decoder:
  decoder <- read.table("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/decoder.txt", header=T,stringsAsFactors=F)

#The count files:
  countFiles <- paste0("/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/",
                      decoder$sample.ID,
                      "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")
# gff
  dexseq.gff <- "/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.flatgff"

# sample table
  sampleTable <- data.frame( row.names = decoder$sample.ID, condition = decoder$group.ID)

#dexseq
  dxd <- DEXSeqDataSetFromHTSeq(
            countFiles,
            sampleData = sampleTable,
            design = ~sample + exon + condition:exon,
            flattenedfile=dexseq.gff)
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd);
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition"); dxres <- DEXSeqResults(dxd);
  dxres
