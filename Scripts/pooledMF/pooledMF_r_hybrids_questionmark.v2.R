#ijob -c1 -p standard -A berglandlab
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

# Question:
## Are the pooled F1s from 2018 D8 basically an even mixture of F1s between A & C?

# Answser: There are two ways to test this.
## The first asks what fraction of polymorphisms in the pooled data exist as informative markers between A & C.
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/old/pooledMF_r_hybrids_questionmark.R` implements the first approach
## The second is to calculate the average frequency of the unique gneotypes among the A & C
### `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/old/pooledMF_r_hybrids_questionmark.R` implements the first approach

### load

### pooled data
  load("/nv/vol186/bergland-lab/alan/totalADRDlongall.Rdata")

  geno[,effRD:=floor(RRD)]
  geno[,effPA:=round(propalt*effRD)/effRD]

  geno.w <- dcast(geno[pond=="D8"], chr + pos ~ Sample, value.var=c("effRD", "effPA"))

### merge
  setkey(genodat.w, chr, pos)
  setkey(geno.w, chr, pos)

  m <- merge(genodat.w, geno.w, all.x=T, all.y=T)

### are pooled D8 samples basically F1s?
  m <- na.omit(m)
  m[,m.hat:=(effPA_D8Male1 * effRD_D8Male1 + effPA_D8Male2 * effRD_D8Male2) / (effRD_D8Male1 + effRD_D8Male2)]
  m[,f.hat:=(effPA_D8PE1 * effRD_D8PE1 + effPA_D8PE2 * effRD_D8PE2) / (effRD_D8PE1 + effRD_D8PE2)]
  m[!is.na(A),A.geno := unlist(sapply(m[!is.na(A)]$A, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]
  m[!is.na(B),B.geno := unlist(sapply(m[!is.na(B)]$B, function(x) c(0,1,2)[which.min(abs(x-c(0,.5,1)))]))]



  m.ag <- m[,list(mu=c(mean(m.hat), mean(f.hat)),
          sd=c(sd(m.hat), sd(f.hat)),
          pool=c("m", "f")),
      list(A.geno, B.geno)]
  m.ag[,exp.fq:=(A.geno+B.geno)/2]


  save(m.ag, file="/nv/vol186/bergland-lab/alan/mf_expectation.Rdata")
