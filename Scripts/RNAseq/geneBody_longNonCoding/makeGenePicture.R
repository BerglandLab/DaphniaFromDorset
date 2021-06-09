BiocManager::install("ballgown")

library(ballgown)
library(ggbio)

gtf <- "/Users/alanbergland/daphnia_ref/Daphnia.aed.0.6.gff"
grl <- gffReadGR(gtf, splitByTranscript = FALSE, identifier = "transcript_id", sep = "; ")


parents <- unlist(grl$Parent)

parents[which(grepl("Daphnia00787", unlist(grl$Parent)))]

autoplot(grl, aes(type = model), gap.geom = "chevron", which="Daphnia00787")
