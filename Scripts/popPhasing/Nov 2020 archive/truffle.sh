#!/usr/bin/env bash
#
##SBATCH -J truffle_pulexOnly # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-03:00:00 # Running time of 5 days
#SBATCH --mem 100G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/truffle_pulexOnly.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/truffle_pulexOnly.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# ijob -c1 -p standard -A berglandlab
### run with: sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/truffle.sh
### sacct -u aob2x -j 11797847
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.10096702_1.err

module load samtools htslib bcftools/1.9

### pulex only samples
  cat /project/berglandlab/Karen/MappingDec2019/CloneInfoFilePulexandObtusa_withmedrd_update20200324 | grep "pulex" | cut -f1 > \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/pulexSamples.list

### make pulex only file
 #bcftools view \
 #-S /scratch/aob2x/daphnia_hwe_sims/popPhase/pulexSamples.list \
 #-Ob \
 #--threads 20 \
 #/scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.vcf.gz > \
 #/scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.pulexOnly.vcf.gz

 bcftools view \
 -Oz \
 --threads 20 \
 /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.pulexOnly.vcf.gz > \
 /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.pulexOnly.vcf.use.gz





### truffle
  twd="/scratch/aob2x/daphnia_hwe_sims/popPhase/truffle"

  ${twd}/truffle \
  --vcf /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.pulexOnly.vcf.use.gz \
  --maf 0.15 \
  --cpu 20  \
  --segments





### library
  library(data.table)
  library(ggplot2)
  library(ggtern)

  ibd <- fread("MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.12chr.LDprune.renameChr.ibd_king.delim")
  ibd <- fread("truffle.ibd")
  ibd[,kinship:=IBD1/4 + IBD2/2]



ggplot(data=ibd, aes(y=PI_HAT, x=IBD0)) + geom_point()

ggplot() +
geom_point(data=ibd,
            aes(x=IBD0, y=IBD1, z=IBD2, size=.85)) +
        coord_tern(expand=T) + limit_tern(T = 1.05, L = 1.05, R = 1.05)

seg <- fread("truffle.segments")


seg.ag <- seg[,list(.N, len=sum(abs(LENGTH))), list(ID1, ID2, TYPE)]

ggplot(seg.ag, aes(x=len, y=N)) + geom_point() + facet_wrap(~TYPE)
