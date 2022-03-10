processQ <- function(k, samp, run, orderby, n=NULL) {
  #k<-8; run=30; samp=samps$clone; orderby="set"

  qi<-as.data.table(q.list[[k]][[run]])
  qi <- as.data.table(qi)
  qi[,gr:=apply(qi, 1, which.max)]
  qi[,samp:=samp]
  qil <- melt(qi, measure.vars=paste("V", 1:k, sep=""))
  setnames(qil, "samp", "clone")
  qil <- merge(qil, samps, by="clone")
  qil[set=="pulex_Dcat", set:="pulex_DCat"]



  if(orderby=="gr") {
    qil.ag <- foreach(gr.i=unique(qil$gr), .combine="rbind")%do%{
      qil.sub <- qil[gr==gr.i][variable==paste("V", gr.i, sep="")]
      qil.sub[,ord:=rank(value, ties.method="first")]
      qil.sub
    }
  } else if(orderby=="set") {
    qil.ag <- foreach(set.i=unique(qil$set), .combine="rbind")%do%{
      #set.i<-"pulex_D8"
      qil.temp <- qil[set==set.i]
      qil.temp.ave <- qil.temp[,list(mu=mean(value)), list(variable)]
      qil.temp <- merge(qil.temp, qil.temp.ave, by="variable")

      qil.temp[,var.rank:=as.numeric(as.factor(rank(mu, ties.method="average")))]

      qil.temp[,wm:=value*var.rank]

      qil.temp.ag <- qil.temp[,list(gr=unique(gr), sum=sum(wm)), list(clone)]

      keycol <-c("gr", "sum")
      setorderv(qil.temp.ag, keycol)
      qil.temp.ag[,ord:=c(1:dim(qil.temp.ag)[1])]
      qil.temp.ag
    }
  } else if(orderby=="set2") {
    qil.ag <- foreach(set.i=unique(qil$set), .combine="rbind")%do%{
      #set.i<-"pulex_D8"

      #set.i<-unique(qil$set)[10]
      qil.temp <- qil[set==set.i]
      qil.temp.ave <- qil.temp[,list(mu=median(value)), list(variable)]
      qil.temp.ave[,rank:=rank(-1*mu, ties.method="first")]
      qil.temp <- merge(qil.temp, qil.temp.ave, by="variable")

      #message("pre")
      qil.temp.ag <- foreach(r=sort(unique(qil.temp$rank)), .combine="merge")%do%{
        #r=1

        qil.temp.rank <- qil.temp[rank==r, c("clone", "value"), with=F]
        qil.temp.rank <- qil.temp.rank[,withinGrRank:=rank(-1*value, ties.method="first")]
        setnames(qil.temp.rank, "withinGrRank", paste("withinGrRank", r, sep="."))
        setkey(qil.temp.rank, "clone")
        qil.temp.rank[,c("clone", paste("withinGrRank", r, sep=".")), with=F]

      }
      #message("post")
      if(is.null(n)) n<- length(names(qil.temp.ag[,-"clone", with=F]))
      setorderv(qil.temp.ag, names(qil.temp.ag[,-"clone", with=F])[1:n])
      qil.temp.ag[,ord:=c(1:dim(qil.temp.ag)[1])]
      qil.temp.ag[,c("clone", "ord"),with=F]
    }

  } else if(orderby=="hclust") {
    qil.ag <- foreach(set.i=unique(qil$set), .combine="rbind")%do%{
      #set.i<-unique(qil$set)[1]

      qil.temp <- qil[set==set.i]
      qil.temp.w <- dcast(qil.temp[,c("variable", "value", "clone"), with=F], clone~variable)

      if(dim(qil.temp.w)[1]>1) {
        d <- dist(as.matrix(qil.temp.w[,-"clone",with=F]))

        hc <- hclust(d, method="centroid")

        image(as.matrix(qil.temp.w[,-"clone",with=F])[hc$order,])

        qil.temp.ag <- data.table(clone=qil.temp.w$clone,
                                  gr=unique(qil.temp$gr),
                                  ord=hc$order)
      } else {
        qil.temp.ag <- data.table(clone=qil.temp.w$clone,
                                  gr=unique(qil.temp$gr),
                                  ord=1)

      }
      qil.temp.ag
    }
  }

  setkey(qil, clone)
  setkey(qil.ag, clone)

  qil <- merge(qil, qil.ag)
  qil[,ordf:=factor(ord, levels=sort(unique(ord)))]
  qil[,set:=factor(set,
                  levels=c("obtusa_DBunk",
                          "pulicaria_Pond22",
                          "pulex_W6",
                          "pulex_W1",
                          "pulex_D10",
                          "pulex_D8",
                          "pulex_DBunk",
                          "pulex_DCat",
                          "pulex_DOil",
                          "pulex_Dramp"))]
  return(qil)
}

plotQ <- function(k.i) {

  ggplot(k.i, aes(x=ord, y=value, fill=as.factor(variable))) +
  geom_bar(position="stack", stat="identity", width = 1) +
  facet_grid(~set, scale="free_x", space="free") +
  geom_text(data=k.i[SC=="A"], aes(x=ord, y=1.1, label="A"), size=5) +
  geom_text(data=k.i[SC=="C"], aes(x=ord, y=1.1, label="C"), size=5) +

  theme(strip.text.x = element_text(angle=90),
        panel.spacing = unit(.2, "lines"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_continuous(expand = c(0,0)) +
  xlab("") +
  guides(fill=FALSE)

}

plotQ.sp <- function(k.i) {
  k.i[SC=="A",col:="A"]
  k.i[SC=="C",col:="C"]
  setorderv(k.i, "col")

  puli.k <- k.i[Species=="pulicaria"][which.max(value)]$variable
  obt.k <- k.i[Species=="obtusa"][which.max(value)]$variable

  puli.plot <- ggplot(data=k.i[variable==puli.k][!is.na(set)], aes(x=set, y=log10(value), color=col)) +
  geom_jitter(width = 0.25) +
  #geom_jitter(data=k.i[variable==puli.k][SC=="A"], aes(x=set, y=log10(value)), color="blue", size=1, width = 0.25) +
  #geom_jitter(data=k.i[variable==puli.k][SC=="C"], aes(x=set, y=log10(value)), color="red", size=1, width = 0.25) +
  coord_flip() +
  xlab("% D. pulicaria ancestry") +
  theme(legend.position="none")

  obt.plot <- ggplot(data=k.i[variable==obt.k][!is.na(set)], aes(x=set, y=log10(value), color=col)) +
  geom_jitter(width = 0.25) +
  #geom_jitter(data=k.i[variable==puli.k][SC=="A"], aes(x=set, y=log10(value)), color="blue", size=1, width = 0.25) +
  #geom_jitter(data=k.i[variable==puli.k][SC=="C"], aes(x=set, y=log10(value)), color="red", size=1, width = 0.25) +
  coord_flip() +
  #labs(title=obt.k) +
  xlab("% D. obtusa ancestry") +

  theme(legend.position="none")


  puli.plot | obt.plot

}

message("processQ / plotQ / plotQ.sp")
