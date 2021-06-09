######################################################################
########## sript for running SHAPEIt to generate phase data ##########
### TODO: 4/10/2019 this runs on one chromosome right now.
### 5/30/2019: functionalized all of it to generate imputed/phased data for all major chromosomes
### I need to fucntionalize it a bit
######################################################################


########################
#### VCF processing ####
########################

    ### compress and tabix unphased VCF file
        bgzip -@10 -c /mnt/spicy_3/Karen/Sp2017/NewMapping/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.vcf  > \
        /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.vcf.gz
        tabix /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.vcf.gz

    ### filter on good SNPs; then retabix
        nohup tabix -h -T /mnt/spicy_3/AlanDaphnia/outputData/snps.delim \
        /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.vcf.gz | \
        bgzip -@10 -c > \
        /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.vcf.gz &

        tabix /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.vcf.gz


    ### make SUB-vcf files for each chromosome
        getChr () {
            echo ${1}
            tabix -h \
            /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.vcf.gz \
            ${1} | bgzip -@10 -c > \
            /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.${1}.vcf.gz
        }
        export -f getChr
        parallel -j1 --gnu getChr ::: $( cut -f1 /mnt/spicy_3/AlanDaphnia/outputData/snps.delim | sort | uniq )


    ### make file of individual names for D8/DBunk
        #grep -E "D8|DBunk" /mnt/spicy_3/AlanDaphnia/outputData/samleids.delim | \
        #awk '{if($3<.05) print $0}' | cut -f1 >  /mnt/spicy_3/AlanDaphnia/vcf/shapeit/D8_DBunk.inds

    ### clean out tri-allelic sites (why are these still here?)
        cleanTri () {
            echo ${1}
            bcftools view --max-alleles 2 --exclude-types indels \
            /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.${1}.vcf.gz | \
            bgzip -@10 -c > \
            /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.${1}.biallelic.vcf.gz
        }
        export -f cleanTri
        parallel -j1 --gnu cleanTri ::: $( cut -f1 /mnt/spicy_3/AlanDaphnia/outputData/snps.delim | sort | uniq )


#######################
#### generate PIRs ####  ### AOB: what is a PIR?
#######################

    ### get sample ids to use: does a bit of filtering on low coverage individuals
      cat /mnt/spicy_3/AlanDaphnia/outputData/samleids.delim | \
      awk '{if($3<.05) print $0}' | cut -f1 >  /mnt/spicy_3/AlanDaphnia/vcf/shapeit/all.inds

    ### function to generate reformatted PIRs input bamlist file
        get_pir_bamlist () {
            sampName=${1} ### sampName=April_2017_D8_101
            chr=${2} ### chr=Scaffold_1863_HRSCAF_2081
            bamDir=/mnt/spicy_3/Karen/20162017FinalMapping/bams
            outFile=/mnt/spicy_3/AlanDaphnia/vcf/shapeit/${chr}.bamlist

            echo ${sampName}"\t"${2}

            echo -e ${sampName}"\t"$( ls ${bamDir}/${sampName}_finalmap.bam )"\t"${chr} | awk '{if(NF==3) print $0}' >> ${outFile}

        }
        export -f get_pir_bamlist
        rm /mnt/spicy_3/AlanDaphnia/vcf/shapeit/*.bamlist
        parallel --gnu -j1 get_pir_bamlist ::: $( cat /mnt/spicy_3/AlanDaphnia/vcf/shapeit/all.inds ) ::: $( cut -f1 /mnt/spicy_3/AlanDaphnia/outputData/snps.delim | sort | uniq )

    ### extract PIRs
        get_pir () {
            chr=${1}
            bamlist_file=/mnt/spicy_3/AlanDaphnia/vcf/shapeit/${chr}.bamlist
            vcf=/mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.${chr}.vcf.gz
            pir=/mnt/spicy_3/AlanDaphnia/vcf/shapeit/${chr}.pir

            echo ${1}"\t"${2}

            /mnt/spicy_3/AlanDaphnia/software/extractPIRs.v1.r68.x86_64/extractPIRs \
            --bam ${bamlist_file} \
            --vcf ${vcf} \
            --out ${pir}
        }
        export -f get_pir

        nohup parallel --gnu -j1 get_pir ::: $( cut -f1 /mnt/spicy_3/AlanDaphnia/outputData/snps.delim | sort | uniq ) &


#####################
#### run shapeit ####
#####################
    runShapeit () {

      echo ${1}

      /mnt/spicy_3/AlanDaphnia/software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit \
      -assemble \
      --states 20 \
      --force \
      --thread 6 \
      --input-pir /mnt/spicy_3/AlanDaphnia/vcf/shapeit/${1}.pir \
      --input-vcf /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.${1}.biallelic.vcf.gz \
      --include-ind /mnt/spicy_3/AlanDaphnia/vcf/shapeit/all.inds \
      --output-log /mnt/spicy_3/AlanDaphnia/vcf/shapeit/${1}.outputlog \
      -O /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.${1}.biallelic.vcf.phased


      /mnt/spicy_3/AlanDaphnia/software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit \
      -convert \
      --input-haps /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.${1}.biallelic.vcf.phased \
      --output-vcf /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.${1}.biallelic.vcf.phased.vcf
    }
    export -f runShapeit

    nohup parallel --gnu -j1 runShapeit ::: $( cut -f1 /mnt/spicy_3/AlanDaphnia/outputData/snps.delim | sort | uniq ) &


    vcf-concat \
    /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.Scaf*bi*phased.vcf \
    | gzip -c > \
    /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.biallelic.phased.vcf.gz

















### defunct

#### prepare files for chromopainter
#### trim up shapeit file to only include variable sites
#    cat /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.Scaffold_1863_HRSCAF_2081.vcf.phased.haps |
#    awk '{
#        nAlleles=0
#        for(i=6; i<=NF; i++) {
#            nAlleles+=$i
#
#        }
#        if(nAlleles>0) print $0
#    }' > /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.Scaffold_1863_HRSCAF_2081.vcf.phased.clean.haps
#
#    ### convert haps file to chromopainter file
#        /mnt/spicy_3/AlanDaphnia/software/fs_4.0.1/impute2chromopainter.pl \
#        /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.Scaffold_1863_HRSCAF_2081.vcf.phased.clean.haps \
#        /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.Scaffold_1863_HRSCAF_2081.vcf.phased.clean
#
#
#### run finestructure
#/mnt/spicy_3/AlanDaphnia/software/fs_4.0.1/fs_linux_glibc2.3 \
#example_cp.cp -n -phasefiles example_cp.phase -recombfiles example_cp.recombfile -idfile example_cp.ids -s1minsnps 5000 -s3iters 10000 -s4iters 10000 -go
#
#
#/mnt/spicy_3/AlanDaphnia/software/fs_4.0.1/fs_linux_glibc2.3 \
#/mnt/spicy_3/AlanDaphnia/chromopainter_output/daps.cp \
#-n \
#-phasefiles /mnt/spicy_3/AlanDaphnia/vcf/shapeit/totalnewmapBfiltsnps25bpindels_snps_filter_pass_lowGQmiss.SNPFilter.Scaffold_1863_HRSCAF_2081.vcf.phased.phase \
#-idfile /mnt/spicy_3/AlanDaphnia/vcf/shapeit/D8_DBunk.inds \
#-s1minsnps 5000 -s3iters 10000 -s4iters 10000 -go
#
