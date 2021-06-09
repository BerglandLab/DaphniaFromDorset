### Run Whatshapp to phase A & C from quartets. Use consensus calls across highly replicated superclones
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/whatshap.prepare_vcf.R`
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/whatshap.run.sh`

### Run RABBIT to reconstruct phase of F1s
  #git clone https://github.com/chaozhi/RABBIT.git

  ~~Apparently does not need phased parents, orininally built for only phenotyped F1s. Now it does all F1s~~

  ### Run RABBIT
  1. `rabbit.runRabbit.sh` runs RABBIT
      (requires `template.m`, runs `rabbit.formatData.R` & `rabbit.parseHaplotypes.py`, and files `chrs.csv`)

  2. `rabbit.plotRabbit.R` makes some diagnostic analysis

  ### convert RABBIT to
  3. `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.convert_output.vitPath.R` makes `pio` object:
    `/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm/all.vitPath.pio.csv`

### RUN lme4qtl gwas with perms
  1. 


  `rabbit.convert_output.R` formats output data
    outputs to `/scratch/aob2x/daphnia_hwe_sims/Rabbit_phase_10cm_onlyPheno`

  3. `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.formatParentHaplos.R` extracts out phased A & C

### Run rQTL
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rQTL.convertData.R` Converts Rabbit output file to rQTL data; also makes separate phenotype data file
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rQTL.run.R` Runs rQTL

### Phenotypic distribution plots & Fraser's v-test
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/vTest.plotPheno.R`

### haplotypes of parents at QTL
  1. `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/rabbit.formatParentHaplos.R`: value=2==ALT; value=1==REF
  2. `DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/parentalHaplotype.plot.R`








`rqtl.R` Runs rQTL on phenotype data (needs to be worked on a bit)













### Phase A & C parents using quartet phasing with Whatshap
`whatshap.prepare_vcf.R` makes the input VCF file
`whatshap.run.sh` runs Whatshap to make output
