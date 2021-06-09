### libraries
  library(data.table)
  library(regioneR)

### load data

  ### precomputed files
    load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/subFiles.Rdata")
    sc[,year:=tstrsplit(clone, "_")[[2]]]

  ### HWE snp-set as defined by `DaphniaPulex20162017Sequencing/AlanAnalysis/HWE_scripts/getTargetSNPs.analysis.workstation.R`
    load(file="/mnt/spicy_3/AlanDaphnia/LD_HWE_slidingWindow/snpSet.Rdata")

  ### pooled seq output as defined by `DaphniaPulex20162017Sequencing/AlanAnalysis/pooledMF/pooled_MF_gather.R`
    load("/mnt/sammas_storage/bergland-lab/alan/mf_fet.Rdata")

    setkey(o, chr, pos)
    setkey(snp.dt, chr, pos)
    o <- merge(o, snp.dt)

  ### raw output of HWE
    load("/mnt/sammas_storage/bergland-lab/alan/hwe_stat.Rdata")



### basic test

  ### merge
    setkey(o, chr, pos)
    setkey(snp.set, chr, pos)

    m <- merge(o[pairs%in%c("D8Male1_D8PE1", "D8Male1_D8PE2", "D8Male2_D8PE1", "D8Male2_D8PE2")], snp.set[py=="D8.2016.2017.2018.2019"], all.x=T, allow.cartesian=TRUE)

  ### test
    fisher.test(table(m$p<0.001, !is.na(m$use)))

### clump-test

  A <- makeGRangesFromDataFrame(as.data.frame(o[pairs%in%c("D8Male1_D8PE1", "D8Male1_D8PE2", "D8Male2_D8PE1", "D8Male2_D8PE2")][use.chr==T][p<1e-3][,list(start=pos), list(chr, end=pos)]),
                           start.field="start", end.field="end")

  B <- makeGRangesFromDataFrame(as.data.frame(snp.set[py=="D8.2016.2017.2018.2019"][,list(chr=unique(chr), start=unique(BP1), end=unique(BP2)), list(block.id)]),
                                start.field="start", end.field="end")


  genome <- getGenomeAndMask(genome=snp.dt[(final.use), list(start=min(pos), end=max(pos)), list(chr)], mask=NA)

  pt <- overlapPermTest(A=A, B=B, genome=genome$genome, ntimes=1000)
  plot(pt)


  ol <- as.data.table(overlapRegions(A=A, B=B,  min.bases=1))
  setnames(ol, "chr", "seqnames")

  p <- ggplot() +
  geom_segment(data=as.data.table(genome$genome),
                aes(x=start, xend=end, y=1, yend=1)) +
  facet_wrap(~seqnames, scales="free_x") +
  geom_jitter(data=as.data.table(A),
                aes(x=start, y=2), size=2, color="red", height=.5) +
  geom_segment(data=as.data.table(B),
                aes(x=start, xend=end, y=3, yend=3), size=2, color="blue")

  p
#### alternative plot

  o[,seqnames:=chr]

  setkey(o, chr, pos)
  setkey(snp.dt, chr, pos)
  o <- merge(o, snp.dt)

  ggplot() +
  geom_segment(data=as.data.table(genome$genome),
                aes(x=start, xend=end, y=1, yend=1)) +
  facet_wrap(~seqnames, scales="free_x") +
  geom_segment(data=as.data.table(B),
                aes(x=start, xend=end, y=8, yend=8), size=2, color="blue") +
  geom_point(data=o[pairs%in%c("D8Male1_D8PE1", "D8Male1_D8PE2", "D8Male2_D8PE1", "D8Male2_D8PE2")][p<1e-4][final.use==T][,list(min.p=min(p), pair=pairs[which.min(p)]), list(seqnames, pos)],
              aes(x=pos, y=-log10(min.p), color=pair))


#### more straight forward test?
  setkey(hwe.stat, chr, pos)
  setkey(o, chr, pos)
  m <- merge(hwe.stat[pond=="D8"][year=="2016.2017.2018.2019"],
             o[pairs%in%c("D8Male1_D8PE1", "D8Male1_D8PE2", "D8Male2_D8PE1", "D8Male2_D8PE2")], all.y=T)


  table(m$p.x<0.05, m$p.y<0.05)




















    ggplot() +
    geom_segment(data=as.data.table(genome$genome),
                  aes(x=start, xend=end, y=1, yend=1)) +
    facet_wrap(~seqnames, scales="free_x") +
    geom_segment(data=as.data.table(B),
                  aes(x=start, xend=end, y=8, yend=8), size=2, color="blue") +
    geom_point(data=o[pairs%in%c("D8Male1_D8PE1", "D8Male1_D8PE2", "D8Male2_D8PE1", "D8Male2_D8PE2")][p<1e-3][final.use==T][,list(min.p=min(p), pair=pairs[which.min(p)]), list(seqnames, pos)],
                aes(x=pos, y=-log10(min.p), color=pair))
