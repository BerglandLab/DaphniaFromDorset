# A. Convert original VCF to BCF and filter on sites that pass quality filtering via one of two files (w/in pulex or across multiple)
  1. for within species `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/filterVCF.within.sh`
  2. for all available daphnids `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/filterVCF.daphnids.sh`
  3. for pulex, but length of daphndis + consensus obtusa / pulicaria `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/filterVCF.hybrid_strategy.sh`

# B. Make jobs file (1 line per chr X clone, downsampled to 1 clone/SC, no lab F1s)
  1. `onePerSC.within.R` this generates 1 file which is the 1 line per chr per SC (passing filter)
  2. `onePerSC.daphnids.R` this generates 1 file which is the 1 line per chr per SC (passing filter)
  3. `onePerSC.hybrid_strategy.R` this generates job file for step C, and also the consensus vcf files for the outgroups

# C. Run Whatshapp to generate read-backed phase-sets
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_SplitVCF.within.sh`
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_SplitVCF.daphnids.sh`
  3. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_SplitVCF.hybrid_strategy.sh`

# D. merge individual vcfs
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/mergeVCF.within.sh`
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/mergeVCF.daphnids.sh`
  3. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/mergeVCF.hybrid_strategy.sh`

# E. Run shapeit on each chromosome separatley
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/shapeit4.within.sh`
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/shapeit4.daphnids.sh`
  3. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/shapeit4.hybrid_strategy.sh`

# F. concatenate the bcfs, and convert to VCF, convert to GDS
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/concatenate_convert.sh`
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/vcf2gds.R`

  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/concatenate_convert.hybrid_strategy.sh`

# G. Analysis (this is quick and dirty analysis using GDS file)
  1. Trees `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.trees.R`

# H. Make FASTA formatted data, one per chromosome
  0. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/makeBed.R` makes additional bed file to filter out based on snps that are not in the filter file
      ```### concatenate bed files
        #sort -k1,1 -k2,2n /scratch/aob2x/daphnia_hwe_sims/popPhase/badSites.bed > /scratch/aob2x/daphnia_hwe_sims/popPhase/badSites.sort.bed
        sed 's/ /\t/g' /scratch/aob2x/daphnia_hwe_sims/popPhase/badSites.bed | bedtools merge -nobuf -i - > \
        /scratch/aob2x/daphnia_hwe_sims/popPhase/badSites.merge.bed

        cat /scratch/aob2x/daphnia_hwe_sims/popPhase/badSites.merge.bed /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.bed | \
        sort -k1,1 -k2,2n | \
        bedtools merge -nobuf -i - > \
        /scratch/aob2x/daphnia_hwe_sims/popPhase/badSites.badRegion.sort.merge.bed

      ```
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/makeFASTA.hybrid_strategy.sh` makes fasta files for each haplotype
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/indexFASTA.hybrid_strategy.sh` indexes the fastas above
  ~~2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/combineFASTA.daphnids.sh` makes 12 fasta, one per chr~~

# I. Build trees on windows
  0. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.makeWindows.trees.R` makes window definition file
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.makeTrees.sh` slices individual FASTA file to make trees and compute Cd. Runs `analysis.ape.trees.R`
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.collectTrees.sh` collects the trees and summarizes

# J. Plots and analysis
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.plotApeTree.R` makes some phylo plots, does some analysis re allele age. Early version
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.polarizeHapolotypes_poolSeq.R` combines the phased haplotypes with the pool-seq to test these hypotheses: (A) A haplotypes should be more common in PE pool & C haplotypes more common in male pool; (B) Dcat_supercloneB should be the same haplotype as male-, ; (C) male+ haplotypes should be more commin in Dbunk





















# F. KING to get basic relationships; is based on observed (not imputed) data; identifies parent offspring pairs
  ### concatenate whatshapp VCF files, rename chromosomes, convert bcf to plink
    1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/bcf2plink.sh`

  ### download KING
    cd /scratch/aob2x/daphnia_hwe_sims/king
    wget http://people.virginia.edu/~wc9c/KING/Linux-king.tar.gz
    tar xzvf Linux-king.tar.gz

  ### run KING
    2. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/king.sh`

# G. duoHMM
  ## fixes phasing errors using pedigree; generates recombination map file


  ### poking about
    module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822
    bcftools view Scaffold_2373_HRSCAF_2879.whatshapp.shapeit.bcf | grep -v "##" | less -S

# duoHMM
  wget https://mathgen.stats.ox.ac.uk/genetics_software/duohmm/duohmm_v0.1.7.tar.gz
  wget https://github.com/molgenis/systemsgenetics/releases/download/1.4.0_20-8.1/GenotypeHarmonizer-1.4.23-dist.tar.gz

  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/duoHMM.sh`

  cp /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.bcf \
  /project/berglandlab/alan.




A few notes

### update bai files (I think they must get copied before the bams which throws the error re: their age)
  #for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/bams/*.bai; do
  #  touch -c ${f}
  #done

### one of the bams is problematic, remove from analysis
  #samtools quickcheck -v /scratch/aob2x/daphnia_hwe_sims/popPhase/bams/*.bam > ~/bad_bams.fofn
  #
  #samtools quickcheck /project/berglandlab/Karen/MappingDec2019/bams/PulexBams/Spring_2017_DBunk_340_finalmap_mdup.bam
  #
  #grep -m1 "#CHROM" /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf | \
  #tr '\t' '\n' | grep "Spring_2017_DBunk_340"

### Karen says that these are the same individual
  #Spring_2016_D8_8.1
  #Spring_2016_D8_8.10

### get 12 big chromosomes
#cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' \
#| sort | uniq | grep -v "chr" | sed 's/^/--chromosome /g' | sed 's/$/ \\/g'
