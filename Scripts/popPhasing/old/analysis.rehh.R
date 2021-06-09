#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R

### libraries
  library(data.table)
  library(rehh)

### read vcf
hh <- data2haplohh(hap_file = "/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.filter.snpsvarpulex.whatshap.vcf",
                 polarize_vcf = FALSE,
                 vcf_reader = "data.table",
                 chr.name="Scaffold_2217_HRSCAF_2652")

#scan <- scan_hh(hh)
save(hh, file="~/hh.Rdata")


scp aob2x@rivanna.hpc.virginia.edu:~/hh.Rdata ~/.

library(rehh)

load("~/hh.Rdata")

res <- calc_ehh(hh,
                mrk = which(hh@positions==5198221),
                include_nhaplo = TRUE)
furcation <- calc_furcation(hh,
                            mrk = which(hh@positions==5198221))
plot(furcation)
haplen <- calc_haplen(furcation)
plot(haplen)
