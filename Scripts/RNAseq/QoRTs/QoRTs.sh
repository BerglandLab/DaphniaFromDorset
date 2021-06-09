#!/usr/bin/env bash
#
#SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-08:00:00 # Running time of 4 days
#SBATCH --mem 60G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=8 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/QoRTs/QoRTs.sh
### sacct -u aob2x -j 20964200
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20964200_8.err

###### SLURM_ARRAY_TASK_ID=1
wd=/scratch/aob2x/daphnia_hwe_sims/

samp=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/samples )
echo $samp


#wget http://hartleys.github.io/QoRTs/QoRTs.jar

java -jar -Xmx60G ~/QoRTs.jar \
QC \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_star_testAligned.sortedByCoord.out.bam \
/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/${samp}
