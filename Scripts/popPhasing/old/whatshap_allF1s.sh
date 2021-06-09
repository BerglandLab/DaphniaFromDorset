#!/usr/bin/env bash
#
#
#SBATCH -J popPhasing_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 4-0:00:00 # Running time of 2 hours
#SBATCH --mem 100G # Memory request of 12 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run with: sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_allF1s.sh
## sacct -u aob2x -j 9922928
## grep "Processing"  /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_whatshapp.9922928_4294967294.err
## cut -f10- /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.phase.vcf | grep "|" | less -S

#ijob -c1 -p standard -A berglandlab
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools
export PATH=$HOME/.local/bin:$PATH

### run whatshap

### update bai files (I think they must get copied before the bams which throws the error re: their age)
  #for f in /scratch/aob2x/daphnia_hwe_sims/popPhase/bams/*.bai; do
  #  touch -c ${f}
  #done

### one of the bams is problematic, remove from analysis
  #samtools quickcheck -v /scratch/aob2x/daphnia_hwe_sims/popPhase/bams/*.bam > ~/bad_bams.fofn
  #
  #samtools quickcheck /project/berglandlab/Karen/MappingDec2019/bams/PulexBams/Spring_2017_DBunk_340_finalmap_mdup.bam
  #
  #grep -m1 "#CHROM" /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf | \
  #tr '\t' '\n' | grep "Spring_2017_DBunk_340"

### Karen says that these are the same individual
  #Spring_2016_D8_8.1
  #Spring_2016_D8_8.10

### get 12 big chromosomes
#cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' \
#| sort | uniq | grep -v "chr" | sed 's/^/--chromosome /g' | sed 's/$/ \\/g'

whatshap \
phase \
-o  /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.phase.vcf \
--chromosome Scaffold_1863_HRSCAF_2081 \
--chromosome Scaffold_1931_HRSCAF_2197 \
--chromosome Scaffold_2158_HRSCAF_2565 \
--chromosome Scaffold_2217_HRSCAF_2652 \
--chromosome Scaffold_2373_HRSCAF_2879 \
--chromosome Scaffold_6786_HRSCAF_7541 \
--chromosome Scaffold_7757_HRSCAF_8726 \
--chromosome Scaffold_9197_HRSCAF_10753 \
--chromosome Scaffold_9198_HRSCAF_10754 \
--chromosome Scaffold_9199_HRSCAF_10755 \
--chromosome Scaffold_9200_HRSCAF_10757 \
--chromosome Scaffold_9201_HRSCAF_10758 \
/scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April17_2018_D8_Male1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April17_2018_D8_Male2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April17_2018_D8_Male3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April17_2018_D8_Male4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April17_2018_D8_Male5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_101_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_103_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_119_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_125_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_130_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_131_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_132_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_133_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_134_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_135_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_136_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_137_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_141_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_142_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_143_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_147_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_150_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_151_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_157_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_165_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_167_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_175_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_179_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_17_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_183_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_18_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_191_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_201_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_202_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_203_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_205_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_209_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_210_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_211_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_212_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_213_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_214_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_215_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_222_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_223_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_227_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_230_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_239_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_248_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_249_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_251_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_254_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_256_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_27_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_285_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_298_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_304_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_305_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_326_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_32_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_338_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_342_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_349_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_355_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_359_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_360_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_373_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_399_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_47_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_48_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_515R_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_55_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_58_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_60_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_65_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_73_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_77_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_79_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_91_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_D8_9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_103_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_108_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_111_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_112_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_119_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_11_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_121_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_128_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_129_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_131_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_132_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_133_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_137_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_139_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_13_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_142_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_145_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_147_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_148_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_149_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_151_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_152_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_153_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_154_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_160_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_162_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_185_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_19_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_203_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_213_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_21_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_230_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_232_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_233_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_237_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_23_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_248_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_253_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_254_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_263_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_266_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_268_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_26_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_27_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_28_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_297_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_302_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_314_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_36_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_38_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_45_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_50_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_52_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_63_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_65_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_6_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_71_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_76_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_7_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_90_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_91_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_96_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DBunk_9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DCat_10_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DCat_2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DCat_3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DCat_4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DCat_5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DLily_1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DMud_78_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DOil_1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DOil_2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DOil_3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DOil_5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_DOil_7_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_Dramp_10_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_Dramp_14_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_Dramp_17_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_Dramp_6_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April_2017_Dramp_9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_D8_MomPE1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_D8_MomPE2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_D8_MomPE3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_D8_MomPE4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_D8_MomPE5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_Male1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_Male2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_Male3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_Male4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_Male5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_MomPE1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_MomPE2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_MomPE3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_MomPE4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DBunk_MomPE5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_Male1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_Male2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_Male3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_Male4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_Male5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_MomPE1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_MomPE2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_MomPE3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_MomPE4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/April29_2018_DCat_MomPE5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/December17_2018_D8_1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_41_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_45_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_46_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_49_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_50_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_54_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_57_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_62_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_63_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_67_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_70_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_74_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_76_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_79_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_81_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_84_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_85_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_A12_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_A13_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_A14_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_A20_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_A21_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_A4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_A6_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_A7_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Fall_2016_D10_A9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Lab_2019_D8_222Male_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Lab_2019_D8_349Male_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE10_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE11_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE12_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE13_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE14_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE15_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE16_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE17_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE18_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE19_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE20_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE21_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE6_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE7_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE8_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_D8_MomPE9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE10_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE11_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE12_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE13_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE14_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE15_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE16_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE17_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE18_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE19_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE20_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE6_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE7_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE8_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DBunk_MomPE9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE10_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE11_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE12_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE13_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE14_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE15_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE16_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE17_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE18_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE19_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE20_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE6_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE7_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE8_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March15_2019_DCat_MomPE9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March_2018_D8_18030_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March_2018_DCat_18004_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_10_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_11_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_12_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_13_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_14_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_15_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_16_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_17_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_18_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_19_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_20_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_21_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_22_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_23_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_24_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_25_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_26_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_27_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_28_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_29_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_30_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_31_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_32_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_33_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_34_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_35_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_36_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_37_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_38_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_39_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_40_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_41_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_43_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_44_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_45_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_46_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_6_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_7_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_8_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_Male1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_Male2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_D8_Male3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_11_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_12_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_13_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_14_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_15_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_16_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_17_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_19_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_20_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_24_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_25_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_27_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_28_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_29_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_30_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_31_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_32_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_33_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_34_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_35_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_36_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_39_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_6_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_7_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_8_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_DBunk_9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_Dcat_1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/March20_2018_Dcat_2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_503_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_515_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_520_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_521_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_532_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_538_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_539_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_542_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_544_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_548_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_612_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_632_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_634_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_663_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_668_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_673_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_726_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_731_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_731SM_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_735_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_752_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_756_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_758_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_770_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_770SM_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_771_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_773_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_773SM_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_776_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_D8_780_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_509_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_514_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_515_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_519_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_520_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_523_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_525_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_535_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_536_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_540_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_547_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_549_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_550_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_554_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_557_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_558_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_561_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_567_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_569_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_575_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_579_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_585_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_589_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_590_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_591_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/May_2017_DBunk_594_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D10_10.1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D10_10.3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D10_10.4_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D10_10.5_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D10_10.6_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.10_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.12_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.14_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.16_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.17_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.18_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.20_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.21_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.23_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.24_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.25_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.26_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.28_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.29_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.2_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.31_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.32_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.35_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.36_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.8_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_D8_8.9_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_W1_1.1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_W6_6.1_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2016_W6_6.3_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_187_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_204_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_224_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_225_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_262_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_282_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_294_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_327_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_329_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_336_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_339_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_340_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_343_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_350_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_362_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_D8_366_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_101_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_106_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_110_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_113_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_116_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_116SM_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_122_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_125_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_143_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_190_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_198_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_209_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_215_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_217_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_222_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_228_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_242_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_252_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_260_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_277_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_281_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_295_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_298_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_315_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_317_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_320_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_321_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_333_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_338_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_347_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_347SM_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_34_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_353_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_359_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_360_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_363_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_367_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_373_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_378_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_380_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_387_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_399_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_64_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_73_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_73SM_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_83_finalmap_mdup.bam \
/scratch/aob2x/daphnia_hwe_sims/popPhase/bams/Spring_2017_DBunk_99_finalmap_mdup.bam
