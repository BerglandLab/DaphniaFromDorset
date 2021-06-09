### get paf
  scp aob2x@rivanna.hpc.virginia.edu:/scratch/aob2x/daphnia_pafs/D84A.daphplx_gasm16ml.taglesspaf ~/.


library(pafr)
library(data.table)
ali <- read_paf("~/D84A.daphplx16asm.tagless.paf", include_tags=FALSE)

unique(ali$tname)[grepl(2158, unique(ali$tname))]

win <- 0
scaf <- unique(as.data.table(ali)[tname=="Scaffold_2158_HRSCAF_2565"][tstart>=(5189407-win)][tend<=(5212730+win)]$qname)

long_ali <- subset(ali, alen > 1e4 & mapq > 40)

dotplot(ali, order_by="provided", label_seqs = TRUE, ordering=list("scaffold_3", c("Scaffold_2158_HRSCAF_2565"))) +
theme_bw() +
ylim(c((5189407-win), (5212730+win))) +
geom_hline(yintercept=5189407) +
geom_hline(yintercept=5212730)

dotplot(long_ali)

unique(as.data.table(long_ali)[tname=="Scaffold_2158_HRSCAF_2565"][tstart>=(5189407-0)][tend<=(5212730+0)]$qname)

plot_synteny(long_ali, q_chrom="FLTH02000094.1", t_chrom="Scaffold_2158_HRSCAF_2565", centre=TRUE) +
    theme_bw()


plot_coverage(long_ali, fill='qname') +
   scale_fill_brewer(palette="Set1")
