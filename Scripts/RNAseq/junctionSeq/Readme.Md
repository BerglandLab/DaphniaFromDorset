### first, run QoRTs on all samples
`DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/QoRTs/QoRTs.sh`

### next, make Decoder file with scaling factors
`DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/QoRTs/makeSizeFactor_decoder.R`

### next, identify novel splice variants with QoRTs
`DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/QoRTs/QoRTs_novelSplice.sh`

### then, run JunctionSeq
`DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/junctionSeq/junctionSeq.sh`
