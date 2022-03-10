scp aob2x@rivanna.hpc.virginia.edu:~/LDlong.Rdata ~/.

library(data.table)
library(ggplot2)

load("~/LDlong.Rdata")


LDlong.ag.ag <- LDlong.ag[pos1!=pos2,
                              list(r2.mean=mean(cor^2), r2.lci=quantile(cor^2, .025), r2.uci=quantile(cor^2, .975)),
                              list(bothSig=(sig.x==T & sig.y==T),
                                   sameChr=(chr.x==chr.y), term=term.x)]



pl <- ggplot(data=LDlong.ag.ag,
       aes(x=sameChr, y=r2.mean, group=interaction(bothSig, sameChr), color=bothSig)) +
geom_point(position=position_dodge(width = .25)) +
geom_errorbar(aes(ymin=r2.lci, ymax=r2.uci), width=.1,
              position=position_dodge(width = .25)) +
facet_wrap(~term) +
theme_bw() +
xlab("Same Chromosome") +
ylab("Linkage disequlibrium (r^2)")

ggsave(pl, file="DaphniaFromDorset/Figures/SFigure10_LinkageQTL_NEEDED/SFigure10.pdf")
