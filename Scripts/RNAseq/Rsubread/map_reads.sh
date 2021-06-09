#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-03:00:00 # Running time of 4 days
#SBATCH --mem 20G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=1-2 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/map_reads.sh
### sacct -u aob2x -j 20407463
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20387599_1.err

module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3

#SLURM_ARRAY_TASK_ID=2
wd=/scratch/aob2x/daphnia_hwe_sims/

samp=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/samples )
echo $samp

Rscript ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/map_reads.R ${samp}
