#!/bin/sh
#
#SBATCH -J run_GEVA
#SBATCH -c 10
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 
#SBATCH --mem 40G
#SBATCH -o ../slurmOut/run_GEVA.%A_%a.out # Standard output
#SBATCH -e ../slurmOut/run_GEVA.%A_%a.out # Standard error
#SBATCH -p standard
#SBATCH --account jcbnunez


########USAGE##########EXAMPLEs
#-#sbatch \
#-#--array=1-$( cat ./Scaffold_names.txt | wc -l) \
#-#run.GEVA.sh \
#-#Scaffold_names.txt \
#-#Killwood_samples.txt 
#-#
#-#sbatch \
#-#--array=1-$( cat ./Scaffold_names.txt | wc -l) \
#-#run.GEVA.sh \
#-#Scaffold_names.txt \
#-#Dorset_samples.txt 
#-#
#-#sbatch \
#-#--array=1-$( cat ./Scaffold_names.txt | wc -l) \
#-#run.GEVA.sh \
#-#Scaffold_names.txt \
#-#England_samples.txt 

########BEGIN##########
### Load Modules
module load tabix
module load bcftools
module load vcftools

### Link programs
GEVA=/home/yey2sn/software/geva
 
### Link inout VCF
INPUT_VCF=/project/berglandlab/alan/phased_daphnia/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf

### Guide file for the QTLs
GUIDE=$1
SAMPLEIDS=$2

### Performance parameters
CPU=$SLURM_CPUS_ON_NODE

### Biological Parameters
REC=1.60e-08 #Lynch et al. -- doi: https://doi.org/10.1101/2020.03.03.974485 
MUT=5.69e-09 #Keith et al. -- doi: https://doi.org/10.1101/gr.191338.115 
NE=862000 #Lynch et al. -- doi: https://doi.org/10.1534/genetics.116.190611 

### Model Files
INITPROBS=/home/yey2sn/software/geva/hmm/hmm_initial_probs.txt
EMITIONS=/home/yey2sn/software/geva/hmm/hmm_emission_probs.txt

###########################################################################
###########################################################################
# Generate Internal Variables
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

CHR=`awk '{ print $1  }' $GUIDE | sed -n ${SLURM_ARRAY_TASK_ID}p`

echo "now processing CHR" ${CHR} 

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################

# VCFtools portion
vcftools --vcf $INPUT_VCF \
--chr ${CHR} \
--maf 0.01 \
--keep $SAMPLEIDS \
--recode --recode-INFO-all \
--out ${CHR}

# add annotation
bcftools query \
-f '%CHROM\t%POS\t%POS\t%CHROM\_%POS\n' \
${CHR}.recode.vcf > ${CHR}.annotation.txt

#Index the annotation
bgzip ${CHR}.annotation.txt
tabix -s1 -b2 -e2 ${CHR}.annotation.txt.gz

# Add annotation
bcftools annotate \
-a ${CHR}.annotation.txt.gz \
-c CHROM,FROM,TO,ID \
--rename-chrs $GUIDE \
${CHR}.recode.vcf > ${CHR}.SNPannotated.forGEVA.vcf 

# Make bin files
$GEVA/geva_v1beta \
--vcf ${CHR}.SNPannotated.forGEVA.vcf  \
--rec $REC \
--out ${CHR}



### Run the TMRCA loop

for i in `cat ${CHR}.marker.txt \ |
sed '1d' |  \
awk '{print $3}'`
do

echo $i

# Run GEVA
$GEVA/geva_v1beta \
--input ${CHR}.bin \
--out RUN_${i} \
--treads $CPU \
--mut $MUT \
--hmm $INITPROBS $EMITIONS \
--Ne $NE \
--position $i 

TMP=`sed -n '7p' RUN_${i}.sites.txt`

if [ -z  "$TMP" ] 
then
	echo  " RUN_${i} is empty "
else
	#Print to common file
	echo -e \
	${CHR} \
	${i} \
	${TMP} \
	$SAMPLEIDS \
	>> $SAMPLEIDS.TMRCA.txt
fi

#Remove Clutter inside the loop
rm RUN_${i}.err
rm RUN_${i}.log
rm RUN_${i}.pairs.txt
rm RUN_${i}.sites.txt

done

#Remove Clutter outside the loop
rm ${CHR}.recode.vcf
rm ${CHR}.annotation.txt.gz
rm ${CHR}.annotation.txt.gz.tbi
rm ${CHR}.log
rm ${CHR}.err
rm ${CHR}.SNPannotated.forGEVA.vcf
rm ${CHR}.bin
rm ${CHR}.sample.txt
rm ${CHR}.marker.txt

