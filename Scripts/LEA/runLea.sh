#!/usr/bin/env bash
#
#SBATCH -J split_and_run # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 5-0:00:00
#SBATCH --mem 5G
#SBATCH -o /scratch/aob2x/lea.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/lea.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# sbatch ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/LEA/runLea.sh
# sacct -j 19319072
# cat  /scratch/aob2x/lea.19319072
module load intel/18.0 intelmpi/18.0 R/3.6.3

wd=/scratch/aob2x/daphnia_hwe_sims/
Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/LEA/runLea.R
