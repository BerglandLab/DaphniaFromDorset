library("topGO")

### load in annotation sheet
  go <- read.xls("DaphniaPulex20162017Sequencing/AlanFigures/SFigure_12_DifferentialExpression/Daphnia_annotation_PANTHER.xls")
  go <- as.data.table(go)


### go analysis
  go.dt <- go[,c("qseqid", "GO"), with=F]
  go.dt <- go.dt[GO!=""]
  go.dt[,gene:=tstrsplit(qseqid, "-")[[1]]]
  go.dt.ag <- go.dt[,list(GO=paste(GO, sep=";")), list(gene)]

  geneGOid <- foreach(i=1:dim(go.dt.ag)[1])%do%{
    strsplit(go.dt.ag[i]$GO, ";")[[1]]
  }
  names(geneGOid) <- go.dt.ag$gene

  allGenes <- dec$qval
  names(allGenes) <- dec$GeneID

  geneNames <- names(geneGOid)
  head(geneNames)

  myInterestingGenes <- dec[(cnA-cnB)%in%c(-2, 2)]$GeneID
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) <- geneNames
  str(geneList)


  GOdata <- new("topGOdata",
                ontology="MF",
                allGenes=geneList,
                annot=annFUN.gene2GO,
                gene2GO=geneGOid)
  graph(GOdata)
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")


  allRes <- GenTable(GOdata, classic = resultFis,
                  orderBy = "weight", ranksOf = "classic", topNodes = 20)

allRes
