#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-03:00:00 # Running time of 4 days
#SBATCH --mem 40G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/makeFeatureCounts.sh
### sacct -u aob2x -j 20407843
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20407843

module load gcc/7.1.0  openmpi/3.1.4 R/4.0.0

#SLURM_ARRAY_TASK_ID=1
wd=/scratch/aob2x/daphnia_hwe_sims/
Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/makeFeatureCounts.R
