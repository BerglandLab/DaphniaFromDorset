#!/usr/bin/env bash
#
#
#SBATCH -J popPhasing_prepareSplitVCF # A single job name for the array
#SBATCH --ntasks=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=12 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 2:00:00 # Running time of 2 hours
#SBATCH --mem 5G # Memory request of 12 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_prepareSplitVCF.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_prepareSplitVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run with: sbatch /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/whatshap_prepareSplitVCF.sh
## sacct -u aob2x -j 9995600

#ijob -c1 -p standard -A berglandlab
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib gparallel/20170822
export PATH=$HOME/.local/bin:$PATH

### first, bgzip and tabix vcf file
  if [ ! -f /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf.gz ]; then
    echo "bgzip original vcf"
    bgzip \
    -c -i -@ 12 \
    /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf > \
    /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf.gz
  fi

  if [ ! -f /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf.gz.gzi ]; then
    echo "tabix bgzipped vcf"
    tabix -p vcf /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf.gz
  fi


### get 12 big chromosomes, and split VCF file
  chr=$( cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/jobId | cut -f2 -d' ' | sort | uniq | grep -v "chr" )

  extractCHR () {
    line=${1}

    echo "... $line ..."

    tabix \
    -p vcf \
    -h \
    /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf.gz \
    $line |
    bgzip -c -@1 > \
    /scratch/aob2x/daphnia_hwe_sims/popPhase/MapDec19PulexandObtusaC_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.${line}.vcf.gz
  }
  export -f extractCHR

  parallel -j 12 extractCHR ::: ${chr}
