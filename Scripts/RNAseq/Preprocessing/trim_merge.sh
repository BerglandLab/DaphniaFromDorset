#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-02:00:00 # Running time of 4 days
#SBATCH --mem 20G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/rna_seq_fastq.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/rna_seq_fastq.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=1-2 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/trim_merge.sh
### sacct -u aob2x -j 20400548
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/rna_seq_fastq.20373409_

module load gcc/9.2.0 bbmap/38.57

#SLURM_ARRAY_TASK_ID=1
wd=/scratch/aob2x/daphnia_hwe_sims/

samp=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/samples )
echo $samp

fn1=$( ls /project/berglandlab/daphnia_rnaseq/usftp21.novogene.com/raw_data/${samp}/${samp}*_1.fq.gz )
fn2=$( ls /project/berglandlab/daphnia_rnaseq/usftp21.novogene.com/raw_data/${samp}/${samp}*_2.fq.gz )
echo $fn1
echo $fn2


bbduk.sh \
in=${fn1} \
in2=${fn2} \
out=/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}_1.trim.fq.gz \
out2=/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}_2.trim.fq.gz \
stats=/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}.trim.stats \
threads=6 \
usejni=t

bbmerge.sh \
in=/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}_1.trim.fq.gz \
in2=/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}_2.trim.fq.gz \
out=/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}.trim.merge.fq.gz \
outu=/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}.trim.r1.fq.gz \
outu2=/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}.trim.r2.fq.gz \
usejni=f
