#!/usr/bin/env bash
#
#SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-02:00:00 # Running time of 4 days
#SBATCH --mem 80G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/junctionSeq/junctionSeq.sh
### sacct -u aob2x -j 20996452
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20996452_

###### SLURM_ARRAY_TASK_ID=1

module load gcc

/home/aob2x/R-3.3.3/bin/Rscript --vanilla /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/junctionSeq/junctionSeq.R
