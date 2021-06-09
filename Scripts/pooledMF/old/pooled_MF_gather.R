### librarioes
  library(data.table)
  library(foreach)

### load data
  fl <- system("ls /scratch/aob2x/daphnia_hwe_sims/mf/mf_fet_*", intern=T)

  o <- foreach(i=fl)%do%{
    print(i)
    fread(i)
  }
  o <- rbindlist(o)

### save
  save(o, file="/nv/vol186/bergland-lab/alan/mf_fet.Rdata")

### plot
  library(data.table)
  library(ggplot2)

  load("/mnt/sammas_storage/bergland-lab/alan/mf_fet.Rdata")

  ggplot(data=o, aes(x=pairs, y=log2(or))) + geom_boxplot() + coord_flip()

  ggplot(data=o[p<.1], aes(x=log2(or), y=-log10(p))) + geom_point() + facet_wrap(~pairs)
