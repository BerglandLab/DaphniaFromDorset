#module load intel/18.0  intelmpi/18.0 gsl/2.4 R/3.6.3 boost/1.60.0; R

#### post process
library(LEA)
library(foreach)
library(data.table)

load(file="/scratch/aob2x/daphnia_hwe_sims/lea/dap.snmf.Rdata")
projectDir = "/scratch/aob2x/daphnia_hwe_sims/lea/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.snmf"

ce = cross.entropy(dap.snmf, K = 2, run = 2)

ce.dt <- foreach(k=1:20, .combine="rbind", .errorhandling="remove")%do%{
	runs <- list.files(paste(projectDir, "/K", k, sep=""))
	foreach(r=runs, .combine="rbind", .errorhandling="remove")%do%{
		print(paste(k, r, sep=" / "))
		#r<-runs[1]
		ce <- fread(paste(projectDir, "/K", k, "/", r, "/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit_r", gsub("run", "", r), ".", k, ".snmfClass", sep=""), header=F, fill=T)[4,2]
		data.table(k=k, run=as.numeric(gsub("run", "", r)),
								ce=as.numeric(gsub("crossEntropy = ", "", ce)))
	}
}

q.list <- foreach(k=1:20, .errorhandling="remove")%do%{
	runs <- list.files(paste(projectDir, "/K", k, sep=""))
	foreach(r=runs, .errorhandling="remove")%do%{
		print(paste(k, r, sep=" / "))
		#r<-runs[1]
		q <- fread(paste(projectDir, "/K", k, "/", r, "/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit_", gsub("run", "r", r), ".", k, ".Q", sep=""), header=F, fill=T)
		q[,k:=k]
		q[,run:=as.numeric(gsub("run", "", r))]
		q
	}
}




### load sample names
	setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/")

### load metadata file
  samps <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
  samps[,set:=paste(Species, population, sep="_")]


	#samps.vcf <- fread("/scratch/aob2x/daphnia_hwe_sims/daps4lea.vcf", skip="#CHROM", nrows=1)
	samps.vcf <- fread("/scratch/aob2x/daphnia_hwe_sims/dsuite/MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf", skip="#CHROM", nrows=1)
	samps.vcf <- names(samps.vcf)[-(1:9)]
	setkey(samps, clone)

	samps <- samps[J(samps.vcf)]
	#samps <- merge(samps.vcf, samps, all.x=T)
	samps[,x:=1:dim(samps)[1]]

	save(ce.dt, q.list, samps, file="/scratch/aob2x/snmf_out.v2.Rdata")
	save(ce.dt, q.list, samps, file="/scratch/aob2x/daphnia_hwe_sims/lea/snmf_out.v2.Rdata")



scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/snmf_out.Rdata ~/.
