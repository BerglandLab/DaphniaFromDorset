#!/usr/bin/env bash
#
#
#SBATCH -J trioPhase_whatshapp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 2:00:00 # Running time of 2 hours
#SBATCH --mem 12G # Memory request of 12 GB
#SBATCH -o /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/daphnia_hwe_sims/slurmOut/trioPhase_whatshapp.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


#ijob -c1 -p standard -A berglandlab
module load gcc/7.1.0 openmpi/3.1.4 python/3.6.8 anaconda/5.2.0-py3.6 samtools
export PATH=$HOME/.local/bin:$PATH

### run whatshap

sed 's/NA/.\/./g' /scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.vcf > \
/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.noNA.vcf

whatshap \
phase \
--ped /scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.ped \
-o /scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.phase.noNA.vcf \
/scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.noNA.vcf


cat /scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.phase.noNA.vcf | grep -v "##" | cut -f1,2,4,5,10- | awk '{
printf $1","$2","$3","$4","
for(i = 5; i <= NF; i++) {
  split($i, sp1, ":")
  printf sp1[1]
  if(i==NF) printf "\n"
  if(i<NF) printf ","
}
}' > /scratch/aob2x/daphnia_hwe_sims/trioPhase/rQTL_quartet.consensus.header.phase.noNA.csv


cp /scratch/aob2x/daphnia_hwe_sims/trioPhase/testTrio.consensus.header.phase.noNA.csv \
/nv/vol186/bergland-lab/alan/.
