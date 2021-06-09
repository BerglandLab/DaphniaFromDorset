library(data.table)

out <- fread("junctionSeq_outallGenes.expression.data.txt.gz")

out[grepl("Daphnia00787", geneID)]
