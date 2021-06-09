First we used a series of commands to prep files to run KING.
```
#!/bin/bash

### Series of commands to prep files for analysis with KING
sed -i 's/Scaffold_1931_HRSCAF_2197/1/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_9198_HRSCAF_10754/2/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_9199_HRSCAF_10755/3/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_9197_HRSCAF_10753/4/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_9200_HRSCAF_10757/5/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_2373_HRSCAF_2879/6/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_7757_HRSCAF_8726/7/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_6786_HRSCAF_7541/8/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_1863_HRSCAF_2081/9/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_2217_HRSCAF_2652/10/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_9201_HRSCAF_10758/11/g' polyDorsetMAF05.vcf
sed -i 's/Scaffold_2158_HRSCAF_2565/12/g' polyDorsetMAF05.vcf


grep "^#" polyDorsetD10WMAF0005.vcf > polyDorsetMAF05.sorted.vcf && grep -v "^#" polyDorsetMAF05.vcf | \
  sort -V -k1,1 -k2,2n >> polyDorsetMAF05.sorted.vcf

grep -v "^Scaff" polyDorsetMAF05.sorted.vcf > polyDorsetMAF05.sortedB.vcf

mv polyDorsetMAF05.sortedB.vcf polyDorsetMAF05.sorted.vcf

module load plink

plink --vcf polyDorsetMAF05.sorted.vcf --out polyDorsetMAF05.sorted --const-fid


grep -v "^#" polyDorsetMAF05.sorted.vcf | grep -v "^Scaff" > polyDorsetMAF05.sorted_tmp

awk '{print $1, $1"_"$2, "0", $2}' polyDorsetMAF05.sorted_tmp > polyDorsetMAF05.sorted.map
```
Then we ran KING using the following slurm script
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -t 6:00:00
#SBATCH --mem=120000
#SBATCH -p largemem
#SBATCH -A berglandlab

#Run program

~/king -b polyDorsetMAF05.sorted.bed --kinship
```
