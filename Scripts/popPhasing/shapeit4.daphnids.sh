#!/usr/bin/env bash
#SBATCH -J popPhasing_shapeit4 # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 02:00:00 # Running time of 2 hours
#SBATCH --mem 120G # Memory request of 120GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_shapeit4.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_shapeit4.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run with: sbatch --array=1-12 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/shapeit4.daphnids.sh
### sacct -u aob2x -j 19104802
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_shapeit4.19104692_1.err
### run with: sbatch --array=9 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/shapeit4.sh
### sacct -u aob2x -j 10267463
#module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822
#
#
## ijob -c1 -p standard -A berglandlab
# SLURM_ARRAY_TASK_ID=1

chr=$( cat /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.daphnids.delim | cut -f1 | sort | uniq | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )

echo $chr

### run shapeit
module load gcc/9.2.0 shapeit4

  shapeit4 \
  --input /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/${chr}.whatshapp.onePerSC.daphnid.bcf \
  --region ${chr} \
  --use-PS 0.0001 \
  --thread 20 \
  --sequencing \
  --log /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/${chr}.daphnid.log \
  --ibd2-length 5 \
  --ibd2-output /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/${chr}.onPerSC.daphnid.IBD2blacklist.txt.gz \
  --output /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/${chr}.whatshapp.onPerSC.daphnid.shapeit.bcf
