#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-02:10:00 # Running time of 4 days
#SBATCH --mem 40G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/coverageDepth/covergeDepth_DB_data.sh
### sacct -u aob2x -j 20865621
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/map_reads.20799222_1.err

module load samtools gcc/9.2.0  openmpi/3.1.6 python/3.7.7 bedtools/2.29.2

#SLURM_ARRAY_TASK_ID=2
wd=/scratch/aob2x/daphnia_hwe_sims/


#cp /project/berglandlab/Karen/RNAseqBams/D86A_merged_sorted.BAM /scratch/aob2x/D86A_merged_sorted.BAM
samtools sort /scratch/aob2x/D86A_merged_sorted.BAM -o /scratch/aob2x/D86A_merged_sorted.sort.BAM -@10
samtools index /scratch/aob2x/D86A_merged_sorted.sort.BAM


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
/scratch/aob2x/D86A_merged_sorted.sort.BAM \
${chr}:${start_win}-${stop_win} > \
/project/berglandlab/alan/bam_slices/D86A_merged_sorted.RNA.small.sort.BAM

samtools index /project/berglandlab/alan/bam_slices/D86A_merged_sorted.RNA.small.sort.BAM







#scp aob2x@rivanna.hpc.virginia.edu:~/*.small.filter.rg.bam* ~/.
#scp aob2x@rivanna.hpc.virginia.edu:~/d8_515_2.small.filter.rg.bam* ~/.
#scp aob2x@rivanna.hpc.virginia.edu:~/*.small.filter.test.bam* ~/.
scp aob2x@rivanna.hpc.virginia.edu:~/*.small.test.bam* ~/.


#chrLen=$( samtools idxstats /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam | grep ${chr} | cut -f2 )
#chrReads=$( samtools idxstats /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam | grep ${chr} | cut -f3 )

#gzip -c /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa > /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa.gz
#python /scratch/aob2x/dest/DEST/mappingPipeline/scripts/PickleRef.py \
#--ref /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa.gz \
#--output /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.pickle

#samtools mpileup \
#-r ${chr}:${start_win}-${stop_win} \
#-B \
#--excl-flags UNMAP \
#-t DP,AD \
#/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam > ~/foo.txt
#less ~/foo.txt
#q
#module load samtools/0.1.20

#samtools mpileup -AB \
#-r ${chr}:${start_win}-${stop_win} \
#/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam > foo.txt
#cat foo.txt | python /scratch/aob2x/daphnia_hwe_sims/allelecount/allelecount.py


#java -ea -Xmx7g -jar /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/coverageDepth/mpileup2sync.jar \
#--fastq-type sanger \
#--min-qual 0 \
#--input ~/foo.txt \
#--output ~/p1_p2.sync

#less ~/p1_p2.sync

#less -S ~/foo.txt
#less ~/foo.txt | awk -v OFS='\t' '{ if ($4>0 && $5 !~ /[^\^][<>]/ && $5 !~ /\+[0-9]+[ACGTNacgtn]+/ && $5 !~ /-[0-9]+[ACGTNacgtn]+/ && $5 !~ /[^\^]\*/) print $1,$2-1,$2,$3,$4,$5,$6}'

#python /scratch/aob2x/dest/DEST_freeze1/mappingPipeline/scripts/Mpileup2Sync.py \
#--mpileup ~/foo.txt \
#--base-quality-threshold 25 \
#--minIndel 1 \
#--coding 1.8 \
#--ref /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.pickle.ref \
#--output ~/output

#cut -f1,2,4 | sed "s/^/$samp\t${chrLen}\t${chrReads}\t/g" > /scratch/aob2x/daphnia_hwe_sims/rnaseq/coverage/${samp}_${gene}.coveragePos.delim
