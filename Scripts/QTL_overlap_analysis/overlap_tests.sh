#!/usr/bin/env bash
#
#SBATCH -J lme4qtl # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0:30:00 # Running time of 1 hours
#SBATCH --mem 12G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run as
# sbatch --array=1-100 ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_overlap_analysis/overlap_tests.sh
# sacct -j 19524846
# cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/lme4qtl.19523871_1.err

module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3

#SLURM_ARRAY_TASK_ID=1
wd="/scratch/aob2x/daphnia_hwe_sims"

Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_overlap_analysis/overlap_tests.R ${SLURM_ARRAY_TASK_ID} AxC
Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/QTL_overlap_analysis/overlap_tests.R ${SLURM_ARRAY_TASK_ID} CxC
