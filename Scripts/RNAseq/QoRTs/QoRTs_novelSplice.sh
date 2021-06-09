#!/usr/bin/env bash
#
#SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-0:02:00 # Running time of 4 days
#SBATCH --mem 20G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/QoRTs/QoRTs_novelSplice.sh
### sacct -u aob2x -j 20988571
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20964200_8.err


#java -Xmx20G -jar ~/QoRTs.jar \
#makeFlatGff \
#/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf \
#/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.flatgff


java -Xmx20G -jar ~/QoRTs.jar \
mergeNovelSplices \
--minCount 10 \
--minSpan 100 \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/ \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/decoder.txt \
/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/qorts_out/

# zcat withNovel.forJunctionSeq.gff.gz | grep "Daphnia00787" | awk '{print $3"\t"$5-$4}'
