#!/usr/bin/env bash
#
##SBATCH -J popPhasing_mergeVCF # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:00:15 # Running time of 4 days
#SBATCH --mem 100G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# ijob -c1 -p standard -A berglandlab
### run with: sbatch --array=1-12 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/mergeVCF.sh
### sacct -u aob2x -j 10259348
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.10096702_1.err

### load modules
  module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822

## get job
  # SLURM_ARRAY_TASK_ID=2
  chr=$( cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' | sort | uniq | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )

# bgzip vcf files
#for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/*.${chr}.phase.vcf; do
#  echo "File -> $f"
#  bgzip \
#  -c \
#  -@ 20 \
#  -i \
#  ${f} > ${f}.gz
#done

# index bgzippped files
#  for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/*.${chr}.phase.vcf.gz; do
#    #f=/scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/April17_2018_D8_Male1.Scaffold_1931_HRSCAF_2197.phase.vcf.gz
#    echo "File -> $f"
#    tabix \
#    -p vcf \
#    -f \
#    ${f}
#  done

  bcftools \
  merge \
  -l /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/${chr}.list \
  -o  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/${chr}.whatshapp.onePerSC.bcf \
  -O b \
  --threads 20

# index
  bcftools index --threads 20 /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/${chr}.whatshapp.onePerSC.bcf
