Filtering the vcf that includes individual genome sequences from 2016,2017, 2018, and 2019, pulex individuals only.
Removing snps within 10 basepairs of indels, hard filtering according to gatk's recommendations, and setting low quality genotype scores to missing.

#First remove all snps within 10 basepairs of indels.

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 3:00:00
#SBATCH --mem=60000
#SBATCH -p standard
#SBATCH -A berglandlab

module load bcftools

cd /scratch/kbb7sh/Daphnia/MappingDecember2019

parameterFile="/scratch/kbb7sh/Daphnia/MappingDecember2019/ChrScaffoldList"

varA=$(sed "${SLURM_ARRAY_TASK_ID}!d" $parameterFile | cut -f 1 )

echo ${varA}

bcftools filter --SnpGap 10 MapJune2020_${varA}.vcf -o MapJune2020_${varA}_filtsnps10bpindels.vcf

```

Next remove all indels, so the file has only SNPs.

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 6:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/MappingDecember2019

echo running gatk

parameterFile="/scratch/kbb7sh/Daphnia/MappingDecember2019/ChrScaffoldList"

varA=$(sed "${SLURM_ARRAY_TASK_ID}!d" $parameterFile | cut -f 1 )

echo ${varA}

java -Xmx4g -jar /scratch/kbb7sh/Daphnia/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R /scratch/kbb7sh/genomefiles/totalHiCwithallbestgapclosed.fa \
	-V /scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_${varA}_filtsnps10bpindels.vcf \
    -selectType SNP \
    -o MapJune2020_${varA}_filtsnps10bpindels_snps.vcf

```

Next hard filter the SNPs based on GATK recommendations for organisms with no reference SNP panel.

First set the filter.

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 4:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/MappingDecember2019

echo running gatk

parameterFile="/scratch/kbb7sh/Daphnia/MappingDecember2019/ChrScaffoldList"

varA=$(sed "${SLURM_ARRAY_TASK_ID}!d" $parameterFile | cut -f 1 )

echo ${varA}

# Run program

java -Xmx4g -jar /scratch/kbb7sh/Daphnia/GenomeAnalysisTK.jar \
	-T VariantFiltration \
	-R /scratch/kbb7sh/genomefiles/totalHiCwithallbestgapclosed.fa \
	-V /scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_${varA}_filtsnps10bpindels_snps.vcf   \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "my_snp_filter" \
    -o MapJune2020_${varA}_filtsnps10bpindels_snps_filter.vcf

```
Then select variants to keep based on filter.
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/MappingDecember2019

echo running gatk

parameterFile="/scratch/kbb7sh/Daphnia/MappingDecember2019/ChrScaffoldList"

varA=$(sed "${SLURM_ARRAY_TASK_ID}!d" $parameterFile | cut -f 1 )

echo ${varA}

# Run program


java -Xmx4g -jar /scratch/kbb7sh/Daphnia/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R /scratch/kbb7sh/genomefiles/totalHiCwithallbestgapclosed.fa \
	-V /scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_${varA}_filtsnps10bpindels_snps_filter.vcf   \
	-ef \
    -o MapJune2020_${varA}_filtsnps10bpindels_snps_filter_pass.vcf
```

Now set low GQ (scores less than 10) to zero.

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/MappingDecember2019

module load bcftools

parameterFile="/scratch/kbb7sh/Daphnia/MappingDecember2019/ChrScaffoldList"

varA=$(sed "${SLURM_ARRAY_TASK_ID}!d" $parameterFile | cut -f 1 )

echo ${varA}

# Run program

bcftools filter -e "FORMAT/GQ<10" -S "." MapJune2020_${varA}_filtsnps10bpindels_snps_filter_pass.vcf | bcftools view -O v -o MapJune2020_${varA}_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf
```
Use gatk check variants to make an index.

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/MappingDecember2019

parameterFile="/scratch/kbb7sh/Daphnia/MappingDecember2019/ChrScaffoldList"

varA=$(sed "${SLURM_ARRAY_TASK_ID}!d" $parameterFile | cut -f 1 )

echo ${varA}

module load gatk

 java -Xmx4g -jar /scratch/kbb7sh/Daphnia/GenomeAnalysisTK.jar \
   -T ValidateVariants \
   -R /scratch/kbb7sh/genomefiles/totalHiCwithallbestgapclosed.fa \
   -V MapJune2020_${varA}_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
   --validationTypeToExclude ALL
```
Next, merge chromosome VCFs into final VCF

```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 2:00:00
#SBATCH --mem=30000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/MappingDecember2019

module load bcftools

bcftools concat \
-o MapJune2020.vcf \
MapJune2020_Scaffold_1863_HRSCAF_2081_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_1931_HRSCAF_2197_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_2158_HRSCAF_2565_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_2217_HRSCAF_2652_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_2373_HRSCAF_2879_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_6786_HRSCAF_7541_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_7757_HRSCAF_8726_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_9197_HRSCAF_10753_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_9198_HRSCAF_10754_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_9199_HRSCAF_10755_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_9200_HRSCAF_10757_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf \
MapJune2020_Scaffold_9201_HRSCAF_10758_filtsnps10bpindels_snps_filter_pass_lowGQmiss.vcf
```
Now annotate using SnpEFF
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -t 4:00:00
#SBATCH --mem=64000
#SBATCH -p largemem
#SBATCH -A berglandlab


cd /scratch/kbb7sh/Daphnia/MappingDecember2019

java -Xmx4g -jar /project/berglandlab/Karen/genomefiles/snpEff/snpEff.jar dpgenome \
/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020.vcf > \
/scratch/kbb7sh/Daphnia/MappingDecember2019/MapJune2020_ann.vcf
```


