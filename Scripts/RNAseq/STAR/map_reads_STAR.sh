#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-04:00:00 # Running time of 4 days
#SBATCH --mem 40G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=1-8 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/STAR/map_reads_STAR.sh
### sacct -u aob2x -j 20592418
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20592384_6.err

module load star/2.7.2b

#SLURM_ARRAY_TASK_ID=6
wd=/scratch/aob2x/daphnia_hwe_sims/

samp=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/samples )
echo $samp

# less /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starLog.final.out
# ls -lh  /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/*_starLog.final.out
# ls -lh  /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}*star*

# STAR \
# --genomeFastaFiles /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.interleaved.fa \
# --runThreadN 5 \
# --runMode genomeGenerate \
# --sjdbGTFfile /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf \
# --genomeDir /project/berglandlab/daphnia_ref/


STAR \
--genomeDir /project/berglandlab/daphnia_ref/ \
--sjdbGTFfile /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf \
--readFilesIn \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}_1.trim.fq.gz \
/scratch/aob2x/daphnia_hwe_sims/rnaseq/trimmed_reads/${samp}_2.trim.fq.gz \
--readFilesCommand gunzip -c \
--sjdbInsertSave All \
--quantMode GeneCounts \
--outReadsUnmapped Fastq \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_star_test \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 38000000000 \
--runThreadN 20 \
--outFilterMatchNmin 0 \
--outSJfilterReads Unique \
--outSJfilterCountUniqueMin 20 1 1 1 \
--alignIntronMax 25000 \
--outFilterMismatchNmax 20 \
--outFilterType BySJout \
--outFilterIntronStrands RemoveInconsistentStrands \
--outWigType bedGraph \
--outWigStrand Unstranded \
--outMultimapperOrder Random \
--sjdbOverhang 100 \
--twopassMode Basic


#--outFilterScoreMinOverLread 0 \
#--outFilterMatchNminOverLread 0 \
