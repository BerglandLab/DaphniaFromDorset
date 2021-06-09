#!/usr/bin/env bash

#SBATCH -J bcf2plink # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 00:10:00 # Running time of 4 days
#SBATCH --mem 20G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/bcf2plink.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/bcf2plink.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run with: sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/bcf2plink.sh
### sacct -u aob2x -j 10127549
#
## ijob -c1 -p standard -A berglandlab
#SLURM_ARRAY_TASK_ID=1

chr=$( cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' | sort | uniq | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )

module load plink/1.90b6.16 bcftools/1.9


### concatenate vcfs
  #ls /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/*.whatshapp.onePerSC.bcf

  ###
  bcftools \
  concat \
  -O b \
  -o /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_1863_HRSCAF_2081.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_1931_HRSCAF_2197.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_2158_HRSCAF_2565.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_2217_HRSCAF_2652.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_2373_HRSCAF_2879.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_6786_HRSCAF_7541.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_7757_HRSCAF_8726.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_9197_HRSCAF_10753.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_9198_HRSCAF_10754.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_9199_HRSCAF_10755.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_9200_HRSCAF_10757.whatshapp.onePerSC.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/Scaffold_9201_HRSCAF_10758.whatshapp.onePerSC.bcf

### rename chrs
  bcftools \
  annotate \
  --rename-chrs /scratch/aob2x/daphnia_hwe_sims/pedigree/chr_rename.delim \
  -O b \
  -o /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.bcf

  #bcftools view /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.bcf | less -S

### rename samples to SC.id
  sed -i 's/,/ /g' /scratch/aob2x/daphnia_hwe_sims/popPhase/representativeSC.delim
  sed '1d' /scratch/aob2x/daphnia_hwe_sims/popPhase/representativeSC.delim  | awk -F' ' '{print $2" "$1}' > /scratch/aob2x/daphnia_hwe_sims/popPhase/representativeSC.rename.delim

  bcftools \
  reheader \
  -s /scratch/aob2x/daphnia_hwe_sims/popPhase/representativeSC.rename.delim \
  -o /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.renameSamp.bcf \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.bcf

### convert to plink
  plink \
  --bcf /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.renameSamp.bcf \
  --make-bed \
  --allow-extra-chr \
  --maf 0.001 \
  --make-founders \
  --indep 100 10 2
  --out /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.renameSamp.ldPrune
