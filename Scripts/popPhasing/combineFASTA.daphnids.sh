#!/usr/bin/env bash
#
##SBATCH -J popPhasing_mergeVCF # A single job name for the array
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-02:00:00 # Running time of 4 days
#SBATCH --mem 10G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# ijob -c1 -p standard -A berglandlab
### run with: sbatch --array=1-12 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/combineFASTA.daphnids.sh
### sacct -u aob2x -j 19110776
### cat /scratch/aob2x/daphnia_hwe_sims/slurmOut/popPhasing_mergeVCF.19110025_10.err

### load modules
  module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools htslib bcftools/1.9 gparallel/20170822

## get job
  # SLURM_ARRAY_TASK_ID=2
  chr=$( cat /scratch/aob2x/daphnia_hwe_sims/popPhase/jobs.id.daphnids.delim | cut -f1 | sort | uniq | grep -v "chr" | awk -v job=${SLURM_ARRAY_TASK_ID} '{if(NR==job) {print $0}}' )

  echo ${chr}

## set up RAM disk
  ## rm /scratch/aob2x/test/*
  #tmpdir="/scratch/aob2x/test"
  #SLURM_JOB_ID=1
  [ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

  tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}


### extract chromosome from each individual
  getChr () {
    chr=${1}
    file=${2} # file="/scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/Spring_2017_DBunk_360.1.fa"
    tmpdir=${3}

    stem=$( echo ${file} | rev | cut -d'/' -f1 | rev )

    echo ${chr} ${stem}


    echo ${chr} > ${tmpdir}/chr.list

    ~/seqtk/seqtk \
    subseq \
    ${file} \
    ${tmpdir}/chr.list | sed "s/${chr}/${stem};${chr}/g" > ${tmpdir}/${chr}_${stem}

  }
  export -f getChr

  parallel getChr ::: ${chr} ::: $( ls -d  /scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/*.fa ) ::: ${tmpdir}

  #cat ${tmpdir}/*.fa > /scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/${chr}.fasta
  awk 'NR>1 && FNR==1{print ""};1' ${tmpdir}/*.fa > /scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/${chr}.fasta

  rm -fr ${tmpdir}
