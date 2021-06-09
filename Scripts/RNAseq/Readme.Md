### 1. Basic quality control
  `DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/rna_seq_fastqc.sh`

  ```
  module load gcc/9.2.0  openmpi/3.1.6 python/3.7.7
  /home/aob2x/.local/bin/multiqc -f .
  tar czv multiqc_data/ multiqc_report.html > ~/daphnia_rnaseq.multiqc.tar.gz
  ```
  scp aob2x@rivanna.hpc.virginia.edu:~/daphnia_rnaseq.multiqc.tar.gz ~/.

### 2. Adapter removal & paired end merging
  `DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/trim_merge.sh`

### 3. Rsubread mapping
  `DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/map_reads.sh`

### 4. merge the PE and "merged" bam files
  `DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/mergebams.sh`

### 5. DEseq2 analysis
  
