#!/usr/bin/env bash
#SBATCH -J ASE_readcounter
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 20G
#SBATCH -t 0-6:00:00
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.%A_%a.err # Standard error

# submit as: sbatch --array=1-8 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/ASE/ASEReadCounter_STAR.sh
# sacct -j 20462313
# cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.20421005_1.err

module load gatk/4.0.0.0 picard samtools gcc/9.2.0 bedtools/2.29.2 vcftools

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf
#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.allvariant.vcf
#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf

### get exonic SNPs
# vcftools \
# --remove-indv pulicaria \
# --remove-indv obtusa \
# --recode \
# --maf 0.00005 \
# --vcf /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf \
# --stdout |
# bedtools \
# intersect \
# -a - \
# -b /project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gff \
# -header -wa -u > \
# /scratch/aob2x/daphnia_hwe_sims/rnaseq/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.exon.pulexOnly.vcf
#
# gatk IndexFeatureFile \
# -F /scratch/aob2x/daphnia_hwe_sims/rnaseq/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.exon.pulexOnly.vcf



#SLURM_ARRAY_TASK_ID=2
wd=/scratch/aob2x/daphnia_hwe_sims/

samp=$( sed "${SLURM_ARRAY_TASK_ID}q;d" ${wd}/DaphniaPulex20162017Sequencing/AlanAnalysis/RNAseq/samples )
echo $samp


###

 #java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
 #-I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.bam \
 #-O /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam \
 #--RGLB ${samp} \
 #--RGPL Illumina \
 #--RGPU ${samp} \
 #--RGSM ${samp}

 #samtools sort \
 #-o /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}.trim.rg.sort.bam \
 #/scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}.trim.rg.bam

gatk ASEReadCounter \
--I /scratch/aob2x/daphnia_hwe_sims/rnaseq/bam/${samp}_starAligned.sortedByCoord.out.rg.bam \
--variant /scratch/aob2x/daphnia_hwe_sims/rnaseq/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.exon.pulexOnly.vcf \
--output /scratch/aob2x/daphnia_hwe_sims/rnaseq/ase/${samp}_rnaseq_asereadcounter.star.delim \
--reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa \
-DF NotDuplicateReadFilter \
