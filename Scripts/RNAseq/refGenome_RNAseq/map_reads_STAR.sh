#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-01:00:00 # Running time of 4 days
#SBATCH --mem 40G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=2-198%1 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/refGenome_RNAseq/map_reads_STAR.sh
### sacct -u aob2x -j 21018001
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.21013208_193.err

module load star/2.7.2b

#SLURM_ARRAY_TASK_ID=6
wd=/scratch/aob2x/daphnia_hwe_sims/

samp=$( ls /project/berglandlab/alan/refGenome_RNAseq/fastq | sed -E 's/_[1-2]{1}\.fq\.gz//g' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
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
/project/berglandlab/alan/refGenome_RNAseq/fastq/${samp}_1.fq.gz \
/project/berglandlab/alan/refGenome_RNAseq/fastq/${samp}_2.fq.gz \
--readFilesCommand gunzip -c \
--sjdbInsertSave All \
--quantMode GeneCounts \
--outReadsUnmapped Fastq \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix /scratch/aob2x/refGenome_RNAseq/bam/${samp}_star \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 38000000000 \
--runThreadN 5 \
--outFilterMatchNmin 0 \
--outSJfilterReads Unique \
--outSJfilterCountUniqueMin 20 1 1 1 \
--alignIntronMax 25000 \
--outFilterMismatchNmax 20 \
--outFilterType BySJout \
--outFilterIntronStrands RemoveInconsistentStrands \
--outWigType bedGraph \
--outWigStrand Stranded \
--outMultimapperOrder Random \
--sjdbOverhang 100

# mv /scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.out.bam \
# /scratch/aob2x/refGenome_RNAseq/bamUse/
#
# mv /scratch/aob2x/refGenome_RNAseq/bam/${samp}_starReadsPerGene.out.tab \
# /scratch/aob2x/refGenome_RNAseq/bamUse/
#
# mv /scratch/aob2x/refGenome_RNAseq/bam/${samp}_starLog.final.out \
# /scratch/aob2x/refGenome_RNAseq/bamUse/
#
# rm -fr /scratch/aob2x/refGenome_RNAseq/bam/${samp}*

#\
#--twopassMode Basic
#test

#--outFilterScoreMinOverLread 0 \
#--outFilterMatchNminOverLread 0 \
