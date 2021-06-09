#!/usr/bin/env bash

#SBATCH -J popPhasing_shapeit4 # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 2-00:00:00 # Running time of 4 days
#SBATCH --mem 100G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_shapeit4.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_shapeit4.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run with: sbatch --array=1-12 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/shapeit4.sh
### sacct -u aob2x -j 10127549
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_shapeit4.10041462_12.err
#module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822
#
#
## ijob -c1 -p standard -A berglandlab
#SLURM_ARRAY_TASK_ID=1
#chr=Scaffold_2158_HRSCAF_2565

chr=$( cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' | sort | uniq | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )

duohmm_wd=/scratch/aob2x/daphnia_hwe_sims/popPhase/duohmm_v0.1.7
gh_wd=/scratch/aob2x/daphnia_hwe_sims/popPhase/GenotypeHarmonizer-1.4.23


module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822 java

bcftools view /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/${chr}.whatshapp.shapeit.bcf > \
/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/${chr}.whatshapp.shapeit.vcf

bgzip -c /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/${chr}.whatshapp.shapeit.vcf > \
/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/${chr}.whatshapp.shapeit.vcf.gz

tabix -p vcf \
/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/${chr}.whatshapp.shapeit.vcf.gz
