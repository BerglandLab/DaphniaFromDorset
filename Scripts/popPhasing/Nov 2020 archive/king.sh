#!/usr/bin/env bash
#
#
#SBATCH -J king # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 6:00:00 # Running time of 15 minutes
#SBATCH --mem 9G # Memory request of 4GGB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/king.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/king.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### run primus: sbatch --array=1-12 /scratch/aob2x/daphnia_hwe_sims/DaphniaPulex20162017Sequencing/AlanAnalysis/popPhasing/king.sh
### sacct -j

kingWD=/scratch/aob2x/daphnia_hwe_sims/king


module load plink/1.90b6.16

### chr=Scaffold_1863_HRSCAF_2081

${kingWD}/king \
-b /scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/12chrs.whatshapp.onePerSC.shapeit.renameChr.renameSamp.ldPrune.bed \
--build \
--prefix /scratch/aob2x/daphnia_hwe_sims/popPhase/kingOut/12chrs.whatshapp.onePerSC.renameChr.renameSamp


plink \
--bfile /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr \
--update-ids /scratch/aob2x/daphnia_hwe_sims/popPhase/kingOut/12chrs.whatshapp.onePerSC.renameChr.updateids.txt \
--make-bed \
--double-id \
--out /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.updateids


plink \
--bfile /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.updateids \
--update-parents /scratch/aob2x/daphnia_hwe_sims/popPhase/kingOut/12chrs.whatshapp.onePerSC.renameChr.updateparents.txt \
--make-bed \
--double-id \
--out /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.updateids.updateparents


plink \
--bfile /scratch/aob2x/daphnia_hwe_sims/popPhase/whatshappOut/12chrs.whatshapp.onePerSC.renameChr.updateids.updateparents \
--recode \
--tab \
--out /scratch/aob2x/daphnia_hwe_sims/popPhase/kingOut/12chrs.whatshapp.onePerSC.renameChr.updateids.updateparents


scp aob2x@rivanna.hpc.virginia.edu:ruffle12chrs.whatshapp.onePerSC.renameChr.updateparents.txt .




library(kinship2)
library(data.table)
ped <- fread("12chrs.whatshapp.onePerSC.renameChr.updateparents.txt")

pedAll <- pedigree(id=c(1:dim(ped)[1]), dadid=ped$V3, momid=ped$V4, sex=rep(4, length(ped$V4)), famid=ped$V1, missid=0)
