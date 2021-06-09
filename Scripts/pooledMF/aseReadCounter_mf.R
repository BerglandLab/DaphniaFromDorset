#ijob -c1 -p standard -A berglandlab
#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.0; R

### libraries
  library(data.table)
  library(SeqArray)
  library(foreach)

### load ASEReadCounter data
  D8Male <- fread("/mnt/sammas_storage/bergland-lab/alan/D8Male.pooledAF.aseReadCounter.allvariant.delim")
  D8PE <- fread("/mnt/sammas_storage/bergland-lab/alan/D8PE.pooledAF.aseReadCounter.allvariant.delim")

### load Karen's data
  load(file="~/mInform.Rdata")

### merge
  setnames(D8Male, c("contig", "position"), c("chr", "pos"))
  setnames(D8PE, c("contig", "position"), c("chr", "pos"))

  setkey(D8Male, chr, pos)
  setkey(m.inform, chr, pos)

  m <- merge(D8Male[,c("chr", "pos", "altCount", "totalCount"), with=F], m.inform[,c("chr", "pos", "effRD_D8Male1", "effRD_D8Male2", "effPA_D8Male1", "effPA_D8Male2", "D8Male.f"), with=F])


### plots
  m[,eff:=(totalCount*160)/(totalCount+160 - 1)]
  plot(eff~I(effRD_D8Male1 + effRD_D8Male2), m)
  abline(0,1)
  m[,freq:=altCount/totalCount]
  plot(I(altCount/totalCount)~D8Male.f, m)
