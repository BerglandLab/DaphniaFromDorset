#!/usr/bin/env bash
#
#SBATCH -J dl_novo # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 8:00:00 ### 0.5 hours
##SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/dl_novo.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/dl_novo.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol8083

# sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/dl_novogene.sh
# sacct -j 12634777
cd /project/berglandlab/alan/

wget -r \
--user=X202SC20051465-Z01_06_04_20_7ML7 \
--password=AEMvQ95U \
ftp://128.120.88.251/
