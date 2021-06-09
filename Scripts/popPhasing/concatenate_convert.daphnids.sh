module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822

ls -d /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/Scaffold*daphnid.shapeit.bcf  > /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/bcfs.daphnid.list

for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/Scaffold*daphnid.shapeit.bcf; do
  #f=/scratch/aob2x/daphnia_hwe_sims/popPhase/tmpFiles/April17_2018_D8_Male1.Scaffold_1931_HRSCAF_2197.phase.vcf.gz
  echo "Index file -> $f"
  bcftools index \
  -f \
  ${f}
done

bcftools \
concat \
-f /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/bcfs.daphnid.list \
-O v \
-l \
-o /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.filter.daphnid.whatshap.vcf


bcftools \
view \
-O b \
-o /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.filter.daphnid.whatshap.bcf \
/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.filter.daphnid.whatshap.vcf

bcftools \
index \
/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.filter.daphnid.whatshap.bcf
