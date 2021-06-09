#!/usr/bin/env bash
#
#SBATCH -J overlap_test # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 1:20:00 # Running time of 1 hours
#SBATCH --mem 10G # Memory request of 8 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/overlap_test.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/overlap_test.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run as
# sbatch --array=1-100 ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/lme4qtl.gwas.sh
# sacct -j 19779660 #AxC
# cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/overlap_test.19779660_1.err

module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3

#SLURM_ARRAY_TASK_ID=1
wd="/scratch/aob2x/daphnia_hwe_sims"

#Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/lme4qtl.gwas.R ${SLURM_ARRAY_TASK_ID} AxC
#Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/lme4qtl.gwas.R ${SLURM_ARRAY_TASK_ID} CxC
Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/rQTL/lme4qtl.gwas.R ${SLURM_ARRAY_TASK_ID} all
