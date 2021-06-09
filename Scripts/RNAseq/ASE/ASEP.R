library(ASEP)
library(data.table)
library(patchwork)
load("~/ase_geno_phase.star.Rdata")

ase.geno.phase[,id:=paste(chr, pos, sep="_")]
ase.simple <- ase.geno.phase[ref_dosage==1][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")][,c("genes", "samp", "id", "xCount", "totalCount", "superclone"), with=F]
setnames(ase.simple, c("genes", "samp", "id", "xCount", "totalCount"), c("gene", "id", "snp", "ref", "total"))
ase.simple[,snp:=as.numeric(as.factor(snp))]


ase.gene.table <- ASE_detection(dat_all = ase.simple[gene%in%paste("Daphnia07379", sep="")],
                                phased=TRUE, varList=NULL, adaptive=TRUE, n_resample=10^3, parallel=FALSE, save_out=FALSE)



dat_M0_phased = as.data.table(phasing(ase.simple[gene=="Daphnia00787"], phased=FALSE, n_condition="one"))

plot_ASE(ase.simple[gene=="Daphnia01430"], phased=T) + theme(legend.position="none")

ggplot(data=ase.simple[gene=="Daphnia07379"], aes(x=id, y=ref/total)) + geom_boxplot()

a1 <- ggplot(data=ase.geno.phase[genes=="Daphnia07379"], aes(x=as.numeric(as.factor(id.x)), y=superclone, fill=as.factor(allele.x))) + geom_tile()
a2 <- ggplot(data=ase.geno.phase[genes=="Daphnia00787"], aes(x=as.numeric(as.factor(id.x)), y=superclone, fill=as.factor(allele.y))) + geom_tile()


ggplot(data=ase.simple[gene%in%paste("Daphnia07379", c(7), sep="")],
              aes(x=id, y=ref/total, fill=clone)) + geom_boxplot() +
              facet_grid(~gene)





load("~/ase_geno_phase.star.Rdata")

ase.simple <- ase.geno.phase[ref_dosage==1][allele.x!=allele.y][superclone=="C"][class%in%c("synonymous_variant", "missense_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")]
ggplot() +
geom_hline( yintercept=mean(ase.simple$xCount/ase.simple$totalCount, na.rm=T)) +
geom_boxplot(data=ase.simple[genes%in%paste("Daphnia0078", c(5:9), sep="")],
              aes(x=clone, y=xCount/totalCount, fill=samp, group=interaction(samp, ref_dosage))) +
facet_grid(~genes, scales="free_x") + theme_bw()

                            ggplot(data=ase.geno.phase[ref_dosage==1][genes%in%paste("Daphnia0078", c(7), sep="")],
                                          aes(x=id, y=xCount/totalCount, fill=superclone, group=clone)) + geom_boxplot() +
                                          facet_grid(~genes, scales="free_x")





ase.simple.ag <- ase.simple[(ref/total)>.05 & (ref/total)<=.95][,
                            list(ref=mean(ref), tot=mean(total), ref_freq=mean(ref)/mean(total), .N),
                            list(gene, id, superclone)]

ase.simple.ag[,p:=pbinom(round(ref), round(tot), .5)]


ggplot(data=ase.simple.ag, aes(x=ref_freq, y=-log10(p))) + geom_point()
ase.simple.ag[ref>0 & ref<.3 & tot>1000]
plot(ref~log10(tot), ase.simple.ag)



a1 / a2
