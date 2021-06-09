# Readme for scripts


# Question:
## Are the pooled F1s from 2018 D8 basically an even mixture of F1s between A & C?

# Answser: There are two ways to test this.
## The first asks what fraction of polymorphisms in the pooled data exist as informative markers between A & C.
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/old/pooledMF_r_hybrids_questionmark.R` implements the first approach

## The second is to calculate the average frequency of the unique gneotypes among the A & C
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/old/pooledMF_r_hybrids_questionmark.v2.R` implements the first approach





## G' test
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/informative_sites.R` makes informative sites info for ASE read counter
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/aseReadCounter_mf.sh` runs ASE read-counter to call allele frequencies
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLSeqR.R` runs Gprime analysis
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/QTLPlot.R` makes composite QLT plot, relies on HARP data





## HWE at QTL regions in F1s
###



find /home/aob2x/ncbi/public/sra -type f -mtime +7 -name '*.cache' -execdir rm -- '{}' \;
