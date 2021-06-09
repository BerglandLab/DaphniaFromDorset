#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:10:00 # Running time of 4 days
#SBATCH --mem 2G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch --array=8-9 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/getBamSlices_Daphnia00787/RNA_bamslices_covergeDepth.sh
### sacct -u aob2x -j 21042935
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.21042935_9.out

module load samtools gcc/9.2.0  openmpi/3.1.6 python/3.7.7 bedtools/2.29.2

#SLURM_ARRAY_TASK_ID=9
wd=/scratch/aob2x/daphnia_hwe_sims/

#samp=$( ls /project/berglandlab/alan/refGenome_RNAseq/fastq | sed -E 's/_[1-2]{1}\.fq\.gz//g' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

#ls /project/berglandlab/alan/refGenome_RNAseq/fastq | sed -E 's/_[1-2]{1}\.fq\.gz//g' > /scratch/aob2x/refJobs
samp=$( cat /scratch/aob2x/refJobs | sed -E 's/_[1-2]{1}\.fq\.gz//g' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

echo ${SLURM_ARRAY_TASK_ID}
echo $samp

if [ ! -f /scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.out.bam.bai ]; then
  samtools index /scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.out.bam
fi

gene=Daphnia00786

chr=$( grep   ${gene} /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff | grep -E "mRNA" | cut -f1 )
start=$( grep ${gene} /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff | grep -E "mRNA" | cut -f4 | tail -n1)
start_win=$( expr $start - 25000 )

gene=Daphnia00789
stop=$( grep  ${gene} /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff | grep -E "mRNA" | cut -f5 | head -n1 )
stop_win=$( expr $stop + 25000 )


#wc -l /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.bed
#cat /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.bed | \
#grep -v "5196681" | grep -v "5201321" | grep -v "5189960" | grep -v "5189615" > /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.keep.bed
#wc -l /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.keep.bed

#samtools view -b \
#/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_star_testAligned.sortedByCoord.out.bam \
#${chr}:${start_win}-${stop_win} |
#bedtools subtract -A -a - -b /project/berglandlab/daphnia_ref/RMoutHiCGMgoodscaff.keep.bed > \
#~/${samp}.small.filter.test.bam
#
#samtools index ~/${samp}.small.filter.test.bam
#

samtools view -b \
/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.out.bam \
${chr}:${start_win}-${stop_win} > \
/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.bam

samtools index /scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.bam


BAM=/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.bam
BAMF1=/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.f1.bam
BAMF2=/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.f2.bam
BAMF=/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.f.bam

BAMR1=/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.r1.bam
BAMR2=/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.r2.bam
BAMR=/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.r.bam

### forward strand
samtools view -bh -f 99 $BAM > $BAMF1
samtools index $BAMF1
samtools view -bh -f 147 $BAM > $BAMF2
samtools index $BAMF2
samtools merge -f $BAMF $BAMF1 $BAMF2
samtools index $BAMF

### reverse strand
samtools view -bh -f 83 $BAM > $BAMR1
samtools index $BAMR1
samtools view -bh -f 163 $BAM > $BAMR2
samtools index $BAMR2
samtools merge -f $BAMR $BAMR1 $BAMR2
samtools index $BAMR



#scp aob2x@rivanna.hpc.virginia.edu:/project/berglandlab/alan/bam_slices/control2_I1.1_D86A_rep1_starAligned.sortedByCoord.small.test.*.ba* ~/.
