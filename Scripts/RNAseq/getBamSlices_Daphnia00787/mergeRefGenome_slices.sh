samtools merge \
-f /scratch/aob2x/refGenome_RNAseq/D86A.forward.bam \
--write-index \
/scratch/aob2x/refGenome_RNAseq/bam/*D86A*_starAligned.sortedByCoord.small.test.f.bam

samtools merge \
-f /scratch/aob2x/refGenome_RNAseq/D86A.reverse.bam \
--write-index \
/scratch/aob2x/refGenome_RNAseq/bam/*D86A*_starAligned.sortedByCoord.small.test.r.bam



samtools merge \
-f /scratch/aob2x/refGenome_RNAseq/C14.forward.bam \
--write-index \
/scratch/aob2x/refGenome_RNAseq/bam/*C14*_starAligned.sortedByCoord.small.test.f.bam

samtools merge \
-f /scratch/aob2x/refGenome_RNAseq/C14.reverse.bam \
--write-index \
/scratch/aob2x/refGenome_RNAseq/bam/*C14*_starAligned.sortedByCoord.small.test.r.bam

scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/refGenome_RNAseq/*forward* ~/bam_slices/
scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/refGenome_RNAseq/*reverse* ~/bam_slices/




BAMR=/scratch/aob2x/refGenome_RNAseq/bam/${samp}_starAligned.sortedByCoord.small.test.r.bam
