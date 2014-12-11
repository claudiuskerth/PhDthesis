## @knitr fragments_mapped_per_ind
setwd("/Users/Claudius/Documents/PhD/THESIS/kks32/Latex/Data_analysis/reference-mapping/data")

fragments_per_ind <- read.delim("fragments_mapped_per_ind",
                            header=T)

par(mfrow=c(1,1), mar=c(4,3,1,3))

mp <- barplot(t(fragments_per_ind[2:5]),
              col=c("red", "yellow", "green", "blue"),
              ylim=c(0,50),
              axes=F
)
points(mp, fragments_per_ind[,length(fragments_per_ind)]/10^5,
       col="violet",
       pch=20,
       type="p"
)
# plot custom x-axis
axis(side=1,
     at=mp,
     labels=fragments_per_ind[,1],
     cex.axis=.6,
     las=2
)
axis(side=4,
     at=seq(0,50,10),
     labels=seq(0,5,1),
     las=1,
     cex.axis=.8
     )
mtext(side=4,
      text="total number of input reads (million)",
      line=2,
      cex=.8
      )
axis(side=2,
     at=seq(0,50,10),
     labels=seq(0,50,10),
     las=1,
     cex.axis=.8
)
mtext(side=2,
      text="number of mapped fragments",
      line=2,
      cex=.8
)

legend("top",
       legend=c("Contig944", "Contig3766", 
                "Contig1776", "Contig213", "total input reads"),
       pch=c(rep(15,4),20),
       col=c("red", "yellow", "green", "blue", "violet"),
       cex=.7,
       pt.cex=1.5,
       bty="n",
       horiz=T,
       text.width=5
)
