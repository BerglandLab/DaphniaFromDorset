#!/usr/bin/env bash
#
##SBATCH -J popPhasing_mergeVCF # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:15git p:00 # Running time of 4 days
#SBATCH --mem 5G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# ijob -c1 -p standard -A berglandlab
### run with: sbatch --array=1-172 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/makeFASTA.hybrid_strategy.sh
### sacct -u aob2x -j 19189296
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.19048068_1.err

### load modules
  module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822 bedtools

## get job
  # SLURM_ARRAY_TASK_ID=2
  #chr=$( cat /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.daphnids.delim | cut -f1 | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )
  samp=$( cat /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.hybrid_strategy.delim | cut -f2 | sort | uniq | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )


###
getHaplos () {
  whichHaplo=${1}
  samp=${2}
  echo ${samp} ${whichHaplo}

  bcftools \
  consensus \
  -f /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa \
  -H ${whichHaplo} \
  -s ${samp} \
  -m /scratch/aob2x/daphnia_hwe_sims/popPhase/badSites.badRegion.sort.merge.bed \
  -o /scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/${samp}.${whichHaplo}.fa \
  /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.bcf

  head /scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/${samp}.${whichHaplo}.fa
}
export -f getHaplos

getHaplos 1 ${samp}
getHaplos 2 ${samp}


#getHaplos 1 pulicaria
#getHaplos 1 obtusa
