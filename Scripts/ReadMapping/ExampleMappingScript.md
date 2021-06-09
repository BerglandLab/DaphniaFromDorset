This is an example of the mapping script that was used to map the Illumina whole genome sequencing data to our European Daphnia pulex reference genome.
```
                ## define some parameters
                sample=${1}
                sampleB=($(echo ${sample} | cut -d"_" -f1-8))
                sampID=($(echo ${sample} | cut -d"_" -f1-7))
                pond=($(echo ${sample} | cut -d"_" -f9-12))
                        
                threads=10
                echo $sample
                echo $sampleB
                inputDir="/scratch/kbb7sh/Daphnia/SingleMoms2018/fastqs"
                interDir="/scratch/kbb7sh/Daphnia/SingleMoms2018/March2018SMIntB"
                outputDir="/scratch/kbb7sh/Daphnia/SingleMoms2018/March2018SMMapB"
                        
                flowcell=($( ls ${inputDir}/*.gz | awk '{split($0,a,"/"); print a[7]}' | awk '{split($0,a,"_"); print a[1]}' | sort | uniq ))
                echo ${!flowcell[@]}            
                                
                        ## map reads
                        for cell in "${flowcell[@]}"; do
                                
                                echo "${cell}"
                                lanes=($( ls ${inputDir}/${cell}*${sampleB}* | awk 'FS="_" {print $2}' | sort | uniq ))
                                        
                                        for lane in "${lanes[@]}"; do
                                                
                                               echo "${lane}"
                                                ### trim out nextera and index seq
                                                java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
                                                        ${inputDir}/${cell}_${lane}_1_${sampleB}.fastq.gz \
                                                        ${inputDir}/${cell}_${lane}_2_${sampleB}.fastq.gz \
                                                        ${interDir}/${cell}_${lane}_1_${sampID}.P_trimm.fastq \
                                                        ${interDir}/${cell}_${lane}_1_${sampID}.U_trimm.fastq \
                                                        ${interDir}/${cell}_${lane}_2_${sampID}.P_trimm.fastq \
                                                        ${interDir}/${cell}_${lane}_2_${sampID}.U_trimm.fastq \
                                                        ILLUMINACLIP:/project/berglandlab/Karen/TrimmomaticAdaptors/NexteraPE-PE.fa:2:30:10:8:true

                                                ### first, merge overlapping reads
                                                ~/pear-0.9.11-linux-x86_64/bin/pear \
                                                -f ${interDir}/${cell}_${lane}_1_${sampID}.P_trimm.fastq \
                                                -r ${interDir}/${cell}_${lane}_2_${sampID}.P_trimm.fastq \
                                                -o ${interDir}/${cell}_${lane}_${sampID} \
                                                -j ${threads}
                
                                                 ### next, map to reference genome
                                                        bwa mem -t ${threads} -K 100000000 -Y \
                                                                -R "@RG\tID:${sampID};${cell};${lane}\tSM:${pond}\tPL:illumina\tPU:${sampID};${cell};${lane}" \
                                                                /scratch/kbb7sh/genomefiles/totalHiCwithallbestgapclosed.fa \
                                                                 ${interDir}/${cell}_${lane}_${sampID}.assembled.fastq | \
                                                        samtools view -L /scratch/kbb7sh/genomefiles/D84Agoodscaffstouse.bed -Suh -q 20 -F 0x100 | \
                                                       samtools sort -@ ${threads} -o ${interDir}/${cell}_${lane}_${pond}.sort.bam
                                                        samtools index ${interDir}/${cell}_${lane}_${pond}.sort.bam

                                                     ## unassembled reads
                                                        bwa mem -t ${threads} -K 100000000 -Y \
                                                                -R "@RG\tID:${sampID};${cell};${lane}\tSM:${pond}\tPL:illumina\tPU:${sampID};${cell};${lane}" \
                                                                /scratch/kbb7sh/genomefiles/totalHiCwithallbestgapclosed.fa \
                                                                ${interDir}/${cell}_${lane}_${sampID}.unassembled.forward.fastq \
                                                                ${interDir}/${cell}_${lane}_${sampID}.unassembled.reverse.fastq | \
                                                        samtools view -L /scratch/kbb7sh/genomefiles/D84Agoodscaffstouse.bed -Suh -q 20 -F 0x100 | \
                                                        samtools sort -@ ${threads} -o ${interDir}/${cell}_${lane}_${pond}.filt.unassembled.sort.bam
                                                        samtools index ${interDir}/${cell}_${lane}_${pond}.filt.unassembled.sort.bam

                                              ## Next, merge assembled and unassembled bam files and mark duplicates
                                                        samtools merge ${interDir}/${cell}_${lane}_${pond}.filt.merged.bam \
                                                                ${interDir}/${cell}_${lane}_${pond}.sort.bam \
                                                                ${interDir}/${cell}_${lane}_${pond}.filt.unassembled.sort.bam
                                                        samtools index ${interDir}/${cell}_${lane}_${pond}.filt.merged.bam
                                                        

                                        done

                        done

                                                ### next, merge bam files to single bam file
                                                samtools merge ${interDir}/${pond}_finalmap.bam ${interDir}/*${pond}.filt.merged.bam
                                                samtools index ${interDir}/${pond}_finalmap.bam 

                                                ### Mark duplicates
                                                java -jar $EBROOTPICARD/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
                                                REMOVE_DUPLICATES=true \
                                                INPUT=${interDir}/${pond}_finalmap.bam \
                                                OUTPUT=${outputDir}/${pond}_finalmap_mdup.bam METRICS_FILE=${interDir}/${pond}_finalmap_mdups.metrics CREATE_INDEX=true
```
The mapping scripts were submitted as arrays using the following example slurm script:
```
#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH -t 8:00:00
#SBATCH --mem=64000
#SBATCH -p standard
#SBATCH -A berglandlab

cd /scratch/kbb7sh/Daphnia/SingleMoms2018

echo mapping

module load gparallel/20170822
module load gcc/7.1.0
module load trimmomatic/0.36
module load bwa
module load samtools
module load picard

#Run program

parameterFile="/scratch/kbb7sh/Daphnia/SingleMoms2018/inputDoerthe"

varA=$(sed "${SLURM_ARRAY_TASK_ID}!d" $parameterFile | cut -f 1 )

bash trialmapping.sh $varA
```
