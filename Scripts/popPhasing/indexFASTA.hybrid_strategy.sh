#!/usr/bin/env bash
#
##SBATCH -J popPhasing_mergeVCF # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:10:00 # Running time of 4 days
#SBATCH --mem 2G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=1-345 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/indexFASTA.hybrid_strategy.sh
### sacct -u aob2x -j 19189468
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.19048068_1.err


# SLURM_ARRAY_TASK_ID=1

### load modules
  module load samtools

file=$( ls -d  /scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/*.fa | sed "${SLURM_ARRAY_TASK_ID}q;d" )

if [ ! -f ${file}.fai ]; then samtools faidx ${file}; fi
