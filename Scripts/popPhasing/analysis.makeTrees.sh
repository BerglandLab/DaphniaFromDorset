#!/usr/bin/env bash
#
##SBATCH -J maketree # A single job name for the array
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1 ### is for multithreading: standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-00:60:00 # Running time of 4 days
#SBATCH --mem 5G # Memory request of 20GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/maketree.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/maketree.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### nJobs=$( cat /scratch/aob2x/daphnia_hwe_sims/popPhase/windows.delim | wc -l ); echo ${nJobs}
### sbatch --array=1-${nJobs} /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.makeTrees.sh
### sacct -u aob2x -j 19791447

### jobs=$( sacct --format ExitCode,JobID%50  -j 19206281 | grep -v "0:0" | cut -f2 -d'_' | sed '1d'| sed '1d'| sed 's/.batch//g' | sed '1d' | tr '\n' ',' | sed 's/ //g' )
### sbatch --array=${jobs} /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.makeTrees.sh
### sacct -j 19210696


### modules
  module load samtools parallel
  module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3

### define parameters
  #chr=Scaffold_2217_HRSCAF_2652
  #start=5173221
  #stop=5223221

  #SLURM_ARRAY_TASK_ID=100

  chr=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/daphnia_hwe_sims/popPhase/windows.delim | cut -f1 -d',' )
  start=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/daphnia_hwe_sims/popPhase/windows.delim | cut -f2 -d',' )
  stop=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/daphnia_hwe_sims/popPhase/windows.delim | cut -f3 -d',' )

  echo ${chr} $start $stop

## set up RAM disk
  ## rm /scratch/aob2x/test/*
  #tmpdir="/scratch/aob2x/test"
  #SLURM_JOB_ID=1; SLURM_ARRAY_TASK_ID=1

  [ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

  tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

### make bed file
  echo ${chr}","${start}","${stop} | tr ',' '\t' > ${tmpdir}/region.bed

### extract chromosome from each individual
  getRegion () {
    chr=${1}
    file=${2} # file="/scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/Spring_2017_DBunk_360.1.fa"
    tmpdir=${3}

    stem=$( echo ${file} | rev | cut -d'/' -f1 | rev )

    echo ${chr} ${stem}

    if [ ! -f ${file}.fai ]; then samtools faidx ${file}; fi

    ~/seqtk/seqtk \
    subseq \
    ${file} \
    ${tmpdir}/region.bed | sed "s/${chr}/${stem};${chr}/g" > ${tmpdir}/${chr}_${stem}


  }
  export -f getRegion

  parallel getRegion ::: ${chr} ::: $( ls -d  /scratch/aob2x/daphnia_hwe_sims/popPhase/FASTA/*.fa ) ::: ${tmpdir}

### combine
  awk 'NR>1 && FNR==1{print ""};1' ${tmpdir}/*.fa > ${tmpdir}/${chr}_${start}_${stop}.fasta

### composition check (for how many Ns)
  ~/seqtk/seqtk \
  comp ${tmpdir}/${chr}_${start}_${stop}.fasta > /scratch/aob2x/daphnia_hwe_sims/popPhase/trees250K/${chr}_${start}_${stop}.info

### make tree
  outgroup=$( grep ">"  ${tmpdir}/${chr}_${start}_${stop}.fasta | grep "obtusa.1.fa" | sed 's/>//g' | sed 's/;/_/g' | sed 's/:/_/g' )

  /home/aob2x/iqtree-1.6.12-Linux/bin/iqtree \
  -s ${tmpdir}/${chr}_${start}_${stop}.fasta \
  -bo 100 \
  -nt 5 \
  -o ${outgroup} \
  -pre ${tmpdir}/${chr}_${start}_${stop}.fasta \
  -m K81 \
  -fast

####
  ### run this R script this does bionj-APE
    #Rscript --vanilla \
    #/scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.ape.trees.R \
    #${tmpdir}/${chr}_${start}_${stop}.fasta \
    #/scratch/aob2x/daphnia_hwe_sims/popPhase/trees250K/${chr}_${start}_${stop}.Rdata

  #cp ${tmpdir}/${chr}_${start}_${stop}.fasta /scratch/aob2x/daphnia_hwe_sims/popPhase/trees/
    Rscript --vanilla \
    /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/analysis.iqtree.trees.R \
    ${tmpdir}/${chr}_${start}_${stop}.fasta \
    /scratch/aob2x/daphnia_hwe_sims/popPhase/trees250K/${chr}_${start}_${stop}.iqtree.Rdata

  rm -fr ${tmpdir}

##

# scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_hwe_sims/popPhase/trees/Scaffold_2217_HRSCAF_2652_5173221_5223221.fasta.treefile ~/.
