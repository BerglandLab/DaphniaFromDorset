#!/usr/bin/env bash
#
##SBATCH -J popPhasing_mergeVCF # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-02:00:00 # Running time of 4 days
#SBATCH --mem 100G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# ijob -c1 -p standard -A berglandlab
### run with: sbatch --array=1-12 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/mergeVCF.daphnids.sh
### sacct -u aob2x -j 19104079
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.19048068_1.err

### load modules
  module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822

## get job
  # SLURM_ARRAY_TASK_ID=2
  chr=$( cat /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.daphnids.delim | cut -f1 | sort | uniq | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )


# bgzip vcf files
for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles.daphnid/*.${chr}.phase.vcf; do
  echo "bgzipping File -> $f"
  bgzip \
  -c \
  -@ 20 \
  -i \
  ${f} > ${f}.gz
done

# index bgzippped files
for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles.daphnid/*.${chr}.phase.vcf.gz; do
  #f=/scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/April17_2018_D8_Male1.Scaffold_1931_HRSCAF_2197.phase.vcf.gz
  echo "indexing File -> $f"
  tabix \
  -p vcf \
  -f \
  ${f}
done

### make file list
  ls -d /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles.daphnid/*.${chr}.phase.vcf.gz > /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles.daphnid/${chr}.list


bcftools \
merge \
-l /scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles.daphnid/${chr}.list \
-o  /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/${chr}.whatshapp.onePerSC.daphnid.bcf \
-O b \
--threads 20

# index
  bcftools index --threads 20 /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/${chr}.whatshapp.onePerSC.daphnid.bcf
