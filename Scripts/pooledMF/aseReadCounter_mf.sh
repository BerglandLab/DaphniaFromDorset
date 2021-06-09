#!/usr/bin/env bash
#SBATCH -J ASE_readcounter
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-4:00:00
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.%A_%a.err # Standard error

# ijob -c1 -p standard -A berglandlab
# submit as: sbatch --array=1-4 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/aseReadCounter_mf.sh
# sacct -j 17864732
# cat /scratch/aob2x/daphnia_hwe_sims/harp_pools/slurmOut/ASE_readcounter.17863582_1.out

module load gatk/4.0.0.0

#cat /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf | awk '{
#  if(substr($1, 0, 1)=="#") {
#    print $0
#  } else {
#    for(i=1; i<=9; i++) printf $i"\t"
#    printf "0/1\t0/1\t0/1\t0/1\t0/1\t0/1\n"
#  }
#}' > /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.allvariant.vcf
#

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.vcf
#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.allvariant.vcf

#gatk IndexFeatureFile -F /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf


###

datadir=/project/berglandlab/Karen/MappingDec2019/MaleFemale2018Pools/bams

if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]; then

  gatk ASEReadCounter \
  --I ${datadir}/HT2LNDSXX_s1_D8PE2.filt.merged.mdup.bam \
  --variant /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf \
  --output /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8PE.pooledAF.aseReadCounter.allvariant.PE2.delim \
  --reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa

elif [[ ${SLURM_ARRAY_TASK_ID} -eq 2 ]]; then

  gatk ASEReadCounter \
  --I ${datadir}/HT2LNDSXX_s1_D8PE1.filt.merged.mdup.bam \
  --variant /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf \
  --output /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8PE.pooledAF.aseReadCounter.allvariant.PE1.delim \
  --reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa


elif [[ ${SLURM_ARRAY_TASK_ID} -eq 3 ]]; then

  gatk ASEReadCounter \
  --I ${datadir}/HT2LNDSXX_s1_D8Male2.filt.merged.mdup.bam \
  --variant /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf \
  --output /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8Male.pooledAF.aseReadCounter.allvariant.Male2.delim \
  --reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa

elif [[ ${SLURM_ARRAY_TASK_ID} -eq 4 ]]; then

  gatk ASEReadCounter \
  --I ${datadir}/HT2LNDSXX_s1_D8Male1.filt.merged.mdup.bam \
  --variant /scratch/aob2x/daphnia_hwe_sims/AC_sites.vcf \
  --output /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8Male.pooledAF.aseReadCounter.allvariant.Male1.delim \
  --reference /project/berglandlab/Karen/MappingDec2019/totalHiCwithallbestgapclosed.fa

fi


SLURM_ARRAY_TASK_ID=1
if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]
then
  echo "foo"
elif [[ ${SLURM_ARRAY_TASK_ID} -eq 2 ]]
then
  echo "bar"
fi


### cp /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8Male.pooledAF.aseReadCounter.allvariant.delim  /nv/vol186/bergland-lab/alan/.
### cp /scratch/aob2x/daphnia_hwe_sims/aseReadCounter/D8PE.pooledAF.aseReadCounter.allvariant.delim  /nv/vol186/bergland-lab/alan/.
