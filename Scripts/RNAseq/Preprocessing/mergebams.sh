#!/usr/bin/env bash
#
##SBATCH -J mergebams # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-01:00:00 # Running time of 4 days
#SBATCH --mem 20G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=1-8 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/mergebams.sh
### sacct -u aob2x -j 20377283
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/rna_seq_fastq.20373409_

module load samtools

#SLURM_ARRAY_TASK_ID=1
wd=/scratch/aob2x/daphnia_hwe_sims/

samp=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/samples )
echo $samp

samtools merge \
${wd}/rnaseq/bam/${samp}.rnaseq.bam \
${wd}/rnaseq/bam/${samp}pe.bam ${wd}/rnaseq/bam/${samp}merge.bam
