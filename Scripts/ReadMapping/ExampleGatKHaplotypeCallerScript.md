This is an example of the scripts that were used to call SNPs for individual clones using gatK's HaplotypeCaller, making individual gvcf files.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH -t 60:00:00
#SBATCH --mem=64000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/MappingDecember2019

echo running gatk 

module load gatk/4.0.0.0

# Run program

parameterFile="/scratch/kbb7sh/Daphnia/MappingDecember2019/oneliterinputs"

varA=$(sed "${SLURM_ARRAY_TASK_ID}!d" $parameterFile | cut -f 1 )

gatk --java-options "-Xmx4g" HaplotypeCaller \
-R /scratch/kbb7sh/genomefiles/totalHiCwithallbestgapclosed.fa \
-I /scratch/kbb7sh/Daphnia/MappingDecember2019/bams/${varA}_finalmap_mdup.bam \
-ERC GVCF \
-O /scratch/kbb7sh/Daphnia/MappingDecember2019/gvcfs/${varA}.g.vcf
```
