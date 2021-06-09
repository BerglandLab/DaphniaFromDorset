#!/usr/bin/env bash
#
#SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-08:00:00 # Running time of 4 days
#SBATCH --mem 1G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/refGenome_RNAseq/getData.sh
### sacct -u aob2x -j 20988852
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20988802


rsync -am \
--include='*.fq.gz' \
--include='*/' \
--exclude='*' \
/nv/vol186/bergland-lab/doerthe/backup_pricey3_19March2020/RNAseq/predation_copper/data/data_A_B_raw/ \
/project/berglandlab/alan/refGenome_RNAseq/fastq
