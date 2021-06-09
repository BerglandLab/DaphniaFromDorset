#!/usr/bin/env bash
#
#SBATCH -J pooledMF # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:15:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/pooledMF.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/pooledMF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load gcc/7.1.0  openmpi/3.1.4 R/3.6.1
Rscript /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/pooled_MF.R \
1000 ${SLURM_ARRAY_TASK_ID}
