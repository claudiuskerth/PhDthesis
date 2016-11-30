# ## ---- set working directory ----
setwd("~/Data_analysis/reference-mapping/data")

## ---- freq dist of genotype calls Big Data ----
# from Big Data standard RAD data set
# assembly with stacks version 0998
# export of all catalog RADtags with their allele depths
# coverage.pl counts genotype calls

rm(list=ls())

geno_calls <- read.delim("export_depth_cov_09042015.tsv", header=T)

# # PyRAD, table 1:
# # my average taxon coverage:
# mean(geno_calls$geno_calls)
# # my % full taxon coverage
# sum(geno_calls$geno_calls == 36)/length(geno_calls$geno_calls) * 100
# # --> these numbers are very similar to the high locus-dropout simulated data set and empircal data set 
# # in the PyRAD paper, table 1
# # number of putative loci:
# length(geno_calls$geno_calls)
# # % singleton loci (i. e. that are only found in one individual):
# sum(geno_calls$geno_calls == 1)/length(geno_calls$geno_calls) * 100



x <- table(geno_calls$geno_calls)
y <- x/sum(x)

# sum(y[1:2])

par(mfrow=c(1,1), bg="cornsilk")

barplot(
        y,
        ylim=c(0, max(y)*1.2),
        ylab="proportion of all loci in the catalog",
        xlab="number of genotype calls",
        main="standard RAD library with SbfI"
)

## ---- freq dist of genotype calls SbfI-XhoI ----
# double digest library with SbfI and XhoI
# merged SE and PE RADtags
# stacks version 09991
# second run with stack
# export of all catalog RADtags with allele depths
# coverage.pl
# from export_sql_allele_depth_20052012.tsv
rm(list=ls())

geno_calls <- read.delim("geno_calls_per_locus", header=T, comment.char="#")
#str(geno_calls)
x <- table(geno_calls$geno_calls)
y <- x/sum(x)

# sum(y[20:38])

par(mfrow=c(1,1), bg="cornsilk")

barplot(
        y,
        ylim=c(0, max(y)*1.2),
        ylab="proportion of all loci in the catalog",
        xlab="number of genotype calls",
        main="double-digest RAD library with SbfI and XhoI"
)

## ---- SbfI frequency dist. per ind. ----

rm(list=ls())

## read in SbfI site positions in uniqued (by individual) PE reads
## see 'position.pl'

## read in SbfI site positions in uniqued (by individual) SE reads

# column number determination on command line:
# $ head -n 1 SE_SbfI_position.out | perl -ne 'print scalar split(","), "\n";'
# $ 90
SE <- read.table("SE_SbfI_position.out", 
                 colClasses = c("character", rep("integer", 89)),
                 sep=",",
                 row.names=1
)

names(SE) <- 1:89
SE <- t(SE)

# column number determination on command line:
# $ head -n 1 PE_SbfI_position.out | perl -ne 'print scalar split(","), "\n";'
# $ 783
# the last field is followed by "\n" not a "," ==> 783 
# $ wc -l PE_SbfI_position.out
# $ 36
PE <- read.table("PE_SbfI_position.out",
                 row.names=1,
                 sep=",",
                 colClasses=c("character", rep("integer", 782))
)

names(PE) <- 1:782
PE <- t(PE)

# Here I am redefining the 'count' function to be able to count
# the occurrence of a single number (n) in a vector (x) OR
# count the occurrences of a vector of numbers (n) in the vector x.
# The function achieves the second mode by letting 'sapply' call
# the 'count' function (i. e. itself) for all numbers in the vector n.

# count <- function(x, n=39){
#         if(length(n) == 1){
#                 return( sum(x==n, na.rm=T) )
#         }else{
#                 return( sapply(X=n, FUN=count, x=x) )
#         }
# }

# For each column in the matrix SE, this returns a vector of length 39
# containing the counts of the numbers 1 to 39 in that column. The 'apply'
# function puts these vectors into a new matrix. The vectors have to be of 
# equal length for the collection into a matrix by the 'apply' function.
# SE.c <- apply(SE, 2, count, n=1:39)
# dim(SE.c)

# Instead of using my new 'count' function, I could have used the 'tabulate'
# function. Note it is important to specify the number of bins (nbins) to get
# a count of 1 to 'nbins' even if no number as high as 'nbins' is found.
SE.c <- apply(SE, 2, tabulate, nbins=39)
# SE.c[,1:5]
row.names(SE.c) <- 1:39
# t(SE.c)[,38:39]

## now, let's plot all SbfI frequency distributions of all indviduals

par(mfrow=c(2,1), mar=c(4,4,3,1))

matplot(SE.c, 
        type="l", 
        col=rgb(1,0,0,.5), 
        lty=1, 
        ylim=c(0,60),
        xlab="position in read (1-based)",
        ylab="unique read count",
        xaxp=c(1,39, 38)
)
mtext("a) single-end reads", side=3, adj=0, line=1)

## now let's do the same for the PE reads

# PE.c <- apply(PE, 2, count, n=1:44)
PE.c <- apply(PE, 2, tabulate, nbins=44)

matplot(PE.c, 
        type="l", 
        col=rgb(0,0,1,.5), 
        lty=1,
        xlab="position in read (1-based)",
        ylab="unique read count",
        xaxp=c(1,44, 43)
)
mtext("b) paired-end reads", side=3, adj=0, line=1)
# There are two individuals with exceptionally high number of SbfI overall.

# class(PE.c)
# dimnames(PE.c)
# colSums(PE.c)
# matplot(PE.c[,c("ery_30-17.fq_2.gz", "ery_30-8.fq_2.gz")], 
#         type="l", 
#         col=rgb(0,0,1,.5), 
#         lty=1,
#         xlab="position in read (1-based)",
#         ylab="unique read count",
#         xaxp=c(1,44, 43)
#         )

# If this is due to fragment-to-fragment religation, then P1-adapter ligation
# should have been inefficient. That should also lead to a reduction in overall
# read number.
# see table 4 and figure 13 in sRAD_Analysis.Rnw
# these two individuals do NOT have an exceptionally low read count

## ---- relative SbfI frequencies ----

PE <- read.table("PE_SbfI_position.out",
                 row.names=1,
                 sep=",",
                 colClasses=c("character", rep("integer", 782))
)

names(PE) <- 1:782
PE <- t(PE)
PE.c <- apply(PE, 2, tabulate, nbins=44)

read.c <- read.delim("read_count", comment.char="#")
#head(read.c)
#head( read.c[order(read.c$file),] )
read.c.sorted <- read.c[order(read.c$file),]
#head(read.c.sorted)
#read.c.sorted$file == dimnames(PE.c)[[2]][order(dimnames(PE.c)[[2]])]
PE.c.sorted <- PE.c[, order(dimnames(PE.c)[[2]])]
#read.c.sorted$file == dimnames(PE.c.sorted)[[2]]
# PE.c[, order(dimnames(PE.c)[[2]])][1:5,1:5]
# PE.c[,grep("ery_30-1\\.|ery_30-10\\.", dimnames(PE.c)[[2]], perl=T)]
# m <- matrix(1:10, 5,2)
# m
# m/c(2,5)
# t(t(m)/c(2,5))
# PE.c.sorted[1:5,1:2]
# read.c.sorted$read_count[1:2]
# t(t(PE.c.sorted[1:5,1:2])/read.c.sorted$read_count[1:2])
PE.c.sorted.rel_freq <- t(t(PE.c.sorted)/read.c.sorted$read_count)

par(mfrow=c(1,1), mar=c(4,4,3,1), bg=gray(0.8))
# class(dimnames(SE.c))
# grep("ery_30-17|ery_30-8", dimnames(SE.c)[[2]])
matplot(PE.c.sorted.rel_freq[,-(grep("ery_30-17|ery_30-8|par_34-3", dimnames(PE.c.sorted.rel_freq)[[2]]))], 
        type="l", 
        col=rgb(0,0,1,0.2), 
        lty=1, 
        ylim=c(0,5e-05),
        xlab="position in read (1-based)",
        ylab="unique read frequency",
        xaxp=c(1,44, 43),
        main="individual frequency distributions of SbfI sites in PE reads"
)
matlines(PE.c.sorted.rel_freq[,(grep("ery_30-17|ery_30-8|par_34-3", dimnames(PE.c.sorted.rel_freq)[[2]]))], 
         type="l", 
         col=c("red", "green", "yellow"), 
         lty=1
)
legend("topright",
       lty=rep(1,3),
       lwd=rep(2,3),
       col=c("red", "green", "yellow"),
       legend=c("ery_30-17", "ery_30-8", "par_34-3")
)

# t(SE.c)[,38:39]

# SE.39 <- apply(SE, 2, count, n=39)
# SE.39
# barplot(sort(SE.39))
# hist(SE.39)
# SE.SbfI.unique.count <- apply(SE, 2, function(x){ sum(!is.na(x)) })
# SE.39.prop <- SE.39/SE.SbfI.unique.count
# barplot(sort(SE.39.prop))
# names(SE.39.prop)
# class(SE.39.prop)
# SE.39.prop.m <- as.matrix(SE.39.prop, ncol=1, nrow=36)
# apply(SE.39.prop.m, 2, sort)
# sort(t(SE.39.prop))
# boxplot(SE)
# 

## ---- fragments_mapped_per_ind ----

fragments_per_ind <- read.delim("fragments_mapped_per_ind",
                            header=T)

#str(fragments_per_ind)

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

## ---- fragments_input_reads_corr ----

#head(fragments_per_ind)
# locus dropout for the 4 contigs by individual
dropout <- apply(fragments_per_ind[,2:5], 1, function(x){sum(x==0)})
reads <- (fragments_per_ind$Retained/10^6) # number of fragments retained after quality fltering and used for analysis

par( mfrow=c(1,1), mar=c(4,4,1,1) )
boxplot(reads ~ dropout, 
     #pch=20, col="red",
     cex.axis=1.2,
     cex.lab=1.2,
     xlim=c(0,4),
     ylim=c(0,5),
     xlab="number of locus-dropout per ind.",
     ylab="input read number per ind. in millions"
)
# simple linear regression
# A GLM would probably be better, but I cannot do that yet.
mod <- lm(reads ~ dropout)
#summary(mod)
R_sq = summary(mod)$r.squared # roughly 16% of the variation in the data can be explained by the linear model
#mod$coeff
# In order to add the regression line to the boxplot,
# I must adjust the intercept to account for the coordinates
# of the boxplot: 0 dropout is at coordinate 1 and so on:
# y = b*x + a
# a = b*1 + ?
# a - b = ?
intercept = mod$coeff["(Intercept)"] - mod$coeff["dropout"]
slope = mod$coeff["dropout"]
abline(intercept, slope, xpd=F, lwd=2)

# add means
mean_reads_per_dropout = aggregate(
        x=list(million_reads=reads), 
        by=list(dropout=dropout), 
        mean)
#mean_reads_per_dropout
points(1:3, mean_reads_per_dropout$million_reads, pch=20)

# add slope, R^squared and the p-value
annot = paste("slope:", round(slope, 3), "(p-value:", "0.0146)")
text(1.5, 4.5, pos=4, cex=1.0,
     #labels=paste(c("slope:", "R^2:", "p-value:"), 
     #            round(c(slope, R_sq, 0.01462), 3)
     labels=annot
)
#round(R_sq, 3)
text(1.5, 4, pos=4, expression(paste(R^2, ": ", 0.163)), cex=1.0)

n = aggregate(x=list(count=rep(1, length(dropout))), 
          by=list(dropout=dropout), 
          FUN=sum)
# class(n)
text(1:3, rep(0, 3), paste(rep("n = ", 3), n[,2]), pos=3)

# fragment count over 4 loci versus input read number
frag_over_loci <- apply(fragments_per_ind[,2:5], 1, sum)
plot(fragments_per_ind$Retained/10^6, frag_over_loci, 
     xlim=c(0,5), 
     cex.lab=1.2,
     cex.axis=1.2,
     ylim=c(0, max(frag_over_loci)+5),
     type="n",
     xlab="total reads input (millions)",
     ylab="sum of mapped fragments (4 loci)"
)
text(fragments_per_ind$Retained/10^6, frag_over_loci, 
     labels=fragments_per_ind[,1], cex=.8)
mod <- lm(frag_over_loci ~ reads)
abline(mod, xpd=F, lwd=2)

## ---- fragNum_per_locus_tab ----

m <- apply(fragments_per_ind[,2:5], 2, mean)
s <- apply(fragments_per_ind[,2:5], 2, sd)
tab <- rbind(m, s/m)
row.names(tab) <- c("mean", "CV")

if (require(xtable)) { cat("")
} else {
        cat("xtable package not found\n")
        cat("installing xtable package now ...\n")
        install.packages("xtable")
        library(xtable)
}

tab <- xtable( t(tab),
        caption="Mean and coefficient of variation of fragment counts per individual for the 4 loci shown in figure \\ref{fragments-mapped-per-ind}.",
        label="mean_sd_fragNum_per_locus",
        display=c("s", "f", "f"),
        digits=c(0, 1, 1),
        align=c("l", "r", "r")
        ) 
print(tab, type="latex", caption.placement="top", booktabs=T)

# kable(
#         t(tab),
#         format = "latex",
#         digits = 1,
#         col.names = c("mean", "CV"),
#         align = c("r", "c", "c"),
#         caption="Mean and coefficient of variation of fragment counts for the 4 loci shown in figure \ref{fragments-mapped-per-ind}.",
#         escape=T,
#         booktabs=T
#         )

        

## ---- ddRAD SbfI frequency dist ----

rm(list=ls())

## read in SbfI site positions in uniqued (by individual) 
## SE reads (see position.pl)
## I deleted the output from the "NotSoTrueTags" files for ery_30-17.

ncol <- max(
        count.fields("ddRAD_SE_SbfI_position.out", sep=",")
)
# ncol
ddRAD_SE_SbfI_pos <- read.csv("ddRAD_SE_SbfI_position.out",
                              sep=",",
                              row.names=1,
                              header=FALSE,
                              fill=TRUE,
                              col.names=paste0("V", seq_len(ncol)),
                              colClasses = c("character", rep("integer", ncol-1)),
)
# str(ddRAD_SE_SbfI_pos)
# row.names(ddRAD_SE_SbfI_pos)

## read in SbfI site positions in uniqued (by individual) 
## PE reads (see position.pl)
ncol <- max(
        count.fields("ddRAD_PE_SbfI_position.out", sep=",")
)
## ncol
ddRAD_PE_SbfI_pos <- read.csv("ddRAD_PE_SbfI_position.out",
                              sep=",",
                              row.names=1,
                              header=FALSE,
                              fill=TRUE,
                              col.names=paste0("V", seq_len(ncol)),
                              colClasses = c("character", rep("integer", ncol-1)),
)

# Instead of using my new 'count' function, I could have used the 'tabulate'
# function. Note it is important to specify the number of bins (nbins) to get
# a count of 1 to 'nbins' even if no number as high as 'nbins' is found.
SE.c <- apply(ddRAD_SE_SbfI_pos, 1, tabulate, nbins=96)
# SE.c[,1:5]
row.names(SE.c) <- 1:96

## now, let's plot all SbfI frequency distributions of all indviduals

## absolute counts

# par(mfrow=c(2,1), mar=c(4,4,3,1))

# matplot(SE.c, 
#         type="l", 
#         col=rgb(1,0,0,.5), 
#         lty=1, 
#         ylim=c(0,60),
#         xlab="position in read (1-based)",
#         ylab="unique read count",
#         xaxp=c(1,96, 95)
# )
# mtext("a) single-end reads", side=3, adj=0, line=1)

## now let's do the same for the PE reads

# PE.c <- apply(PE, 2, count, n=1:44)
PE.c <- apply(ddRAD_PE_SbfI_pos, 1, tabulate, nbins=100)

# matplot(PE.c, 
#         type="l", 
#         col=rgb(0,0,1,.5), 
#         lty=1,
#         xlab="position in read (1-based)",
#         ylab="unique read count",
#         xaxp=c(1,100, 99)
# )
# mtext("b) paired-end reads", side=3, adj=0, line=1)

#
## realtive counts 
#
read.c <- read.delim("record_count", comment.char="#", header=F, sep=" ")
# head(read.c)
names(read.c) <- c("sample", "count")
# head( read.c[order(read.c$sample),] )
read.c.sorted <- read.c[order(read.c$sample),]
#head(read.c.sorted)
#read.c.sorted$sample == dimnames(SE.c)[[2]][order(dimnames(SE.c)[[2]])]
SE.c.sorted <- SE.c[, order(dimnames(SE.c)[[2]])]
# read.c.sorted$sample == dimnames(SE.c.sorted)[[2]]
# diving a matrix by a vector in R is row-wise
# m <- matrix(1:10, 5,2)
# m
# m/c(2,5)
# so I have to turn columns into rows here in order to divide the first column
# by the first element in the vector, the second column by the second element in the vector,
# and so on
# t(t(m)/c(2,5))
# now I can divide the SbfI count for each individual by the number of reads for 
# that individual, thus getting relative counts
SE.c.sorted.rel_freq <- t(t(SE.c.sorted)/read.c.sorted$count)

## now get relative frequencies for PE reads
## note, sample names now end with "fq_2" instead of "fq_1"
## I changed all sample names to end with "fq_1"
PE.c.sorted <- PE.c[, order(dimnames(PE.c)[[2]])]
# read.c.sorted$sample == dimnames(PE.c.sorted)[[2]]
# data.frame(read.c.sorted$sample, dimnames(PE.c.sorted)[[2]])
# dimnames(PE.c.sorted)[[2]]
PE.c.sorted.rel_freq <- t(t(PE.c.sorted)/read.c.sorted$count)

## matplots of relative SbfI frequencies

par(mfrow=c(2,1), mar=c(4,4,3,1))

matplot(SE.c.sorted.rel_freq, 
        type="l", 
        col=rgb(1,0,0,.5), 
        lty=1, 
        #         ylim=c(0, 6e-5),
        xlab="position in read (1-based)",
        ylab="unique read frequency",
        xaxp=c(1,96, 19)
)
mtext("a) single-end reads", side=3, adj=0, line=1)

## now let's do the same for the PE reads

matplot(PE.c.sorted.rel_freq, 
        type="l", 
        col=rgb(0,0,1,.5), 
        lty=1,
        xlab="position in read (1-based)",
        ylab="unique read frequency",
        xaxp=c(1,100, 33)
)
mtext("b) paired-end reads", side=3, adj=0, line=1)

## ---- ddRAD XhoI frequency dist ----

rm(list=ls())

## read in XhoI site positions in uniqued (by individual) 
## SE reads (see position.pl)
## I deleted the output from the "NotSoTrueTags" files for ery_30-17.

ncol <- max(
        count.fields("ddRAD_SE_XhoI_position.out", sep=",")
)
# ncol
ddRAD_SE_XhoI_pos <- read.csv("ddRAD_SE_XhoI_position.out",
                              sep=",",
                              row.names=1,
                              header=FALSE,
                              fill=TRUE,
                              col.names=paste0("V", seq_len(ncol)),
                              colClasses = c("character", rep("integer", ncol-1)),
)


## read in XhoI site positions in uniqued (by individual) 
## PE reads (see position.pl)
ncol <- max(
        count.fields("ddRAD_PE_XhoI_position.out", sep=",")
)
## ncol
ddRAD_PE_XhoI_pos <- read.csv("ddRAD_PE_XhoI_position.out",
                              sep=",",
                              row.names=1,
                              header=FALSE,
                              fill=TRUE,
                              col.names=paste0("V", seq_len(ncol)),
                              colClasses = c("character", rep("integer", ncol-1)),
)

# Instead of using my new 'count' function, I could have used the 'tabulate'
# function. Note it is important to specify the number of bins (nbins) to get
# a count of 1 to 'nbins' even if no number as high as 'nbins' is found.
SE.c <- apply(ddRAD_SE_XhoI_pos, 1, tabulate, nbins=96)
# SE.c[,1:5]
row.names(SE.c) <- 1:96

## now, let's plot all XhoI frequency distributions of all indviduals

## absolute counts

# par(mfrow=c(2,1), mar=c(4,4,3,1))

# matplot(SE.c, 
#         type="l", 
#         col=rgb(1,0,0,.5), 
#         lty=1, 
#         ylim=c(0,60),
#         xlab="position in read (1-based)",
#         ylab="unique read count",
#         xaxp=c(1,96, 95)
# )
# mtext("a) single-end reads", side=3, adj=0, line=1)

## now let's do the same for the PE reads

# PE.c <- apply(PE, 2, count, n=1:44)
PE.c <- apply(ddRAD_PE_XhoI_pos, 1, tabulate, nbins=100)

# matplot(PE.c, 
#         type="l", 
#         col=rgb(0,0,1,.5), 
#         lty=1,
#         xlab="position in read (1-based)",
#         ylab="unique read count",
#         xaxp=c(1,100, 99)
# )
# mtext("b) paired-end reads", side=3, adj=0, line=1)

#
## realtive counts 
#
read.c <- read.delim("record_count", comment.char="#", header=F, sep=" ")
# head(read.c)
names(read.c) <- c("sample", "count")
# head( read.c[order(read.c$sample),] )
read.c.sorted <- read.c[order(read.c$sample),]
#head(read.c.sorted)
#read.c.sorted$sample == dimnames(SE.c)[[2]][order(dimnames(SE.c)[[2]])]
SE.c.sorted <- SE.c[, order(dimnames(SE.c)[[2]])]
# read.c.sorted$sample == dimnames(SE.c.sorted)[[2]]
# diving a matrix by a vector in R is row-wise
# m <- matrix(1:10, 5,2)
# m
# m/c(2,5)
# so I have to turn columns into rows here in order to divide the first column
# by the first element in the vector, the second column by the second element in the vector,
# and so on
# t(t(m)/c(2,5))
# now I can divide the XhoI count for each individual by the number of reads for 
# that individual, thus getting relative counts
SE.c.sorted.rel_freq <- t(t(SE.c.sorted)/read.c.sorted$count)

## now get relative frequencies for PE reads
## note, sample names now end with "fq_2" instead of "fq_1"
## I changed all sample names to end with "fq_1"
PE.c.sorted <- PE.c[, order(dimnames(PE.c)[[2]])]
# read.c.sorted$sample == dimnames(PE.c.sorted)[[2]]
# data.frame(read.c.sorted$sample, dimnames(PE.c.sorted)[[2]])
# dimnames(PE.c.sorted)[[2]]
PE.c.sorted.rel_freq <- t(t(PE.c.sorted)/read.c.sorted$count)

## matplots of relative XhoI frequencies

par(mfrow=c(2,1), mar=c(4,4,3,1))

matplot(SE.c.sorted.rel_freq, 
        type="l", 
        col=rgb(1,0,0,.5), 
        lty=1, 
        #         ylim=c(0, 6e-5),
        xlab="position in read (1-based)",
        ylab="unique read frequency",
        xaxp=c(1,96, 19)
)
mtext("a) single-end reads", side=3, adj=0, line=1)

## now let's do the same for the PE reads

matplot(PE.c.sorted.rel_freq, 
        type="l", 
        col=rgb(0,0,1,.5), 
        lty=1,
        xlab="position in read (1-based)",
        ylab="unique read frequency",
        xaxp=c(1,100, 33)
)
mtext("b) paired-end reads", side=3, adj=0, line=1)

## ---- cluster size by XhoI position ----

rm(list=ls())

ncol <- max(
        count.fields("ddRAD_SE_XhoI_position.out", sep=",")
)
## ncol
ddRAD_SE_XhoI_pos <- read.csv("ddRAD_SE_XhoI_position.out",
                              sep=",",
                              row.names=1,
                              header=FALSE,
                              fill=TRUE,
                              col.names=paste0("V", seq_len(ncol)),
                              colClasses = c("character", rep("integer", ncol-1)),
)

# Instead of using my new 'count' function, I can use the 'tabulate'
# function. Note it is important to specify the number of bins (nbins) to get
# a count of 1 to 'nbins' even if no number as high as 'nbins' is found.
SE.c <- apply(ddRAD_SE_XhoI_pos, 1, tabulate, nbins=96)
# SE.c[91:96,1:5]
# str(SE.c)
SE.c <- SE.c[7:91,]
SE.c.total <- rowSums(SE.c)
# length(7:91)
# length(SE.c.total)

## plot total XhoI site frequency spectrum
par(mar=c(5,4,1,4))

plot(7:91, SE.c.total,
     type="b",
     pch=20,
     col=rgb(1,0,0, 0.5),
     xaxp=c(7, 91, 4),
     xlab="position in read (1-based)",
     ylab="",
     yaxt="n",
     bty ="n"
)
axis(2, col="red")
mtext("unique read count", side=2, col="red", line=3)


ncol <- max(
        count.fields("all_ind_XhoI_in_SE_cl_size_by_pos.csv", sep=",")
)
# ncol
cl <- read.csv("all_ind_XhoI_in_SE_cl_size_by_pos.csv",
               sep=",",
               row.names=1,
               header=FALSE,
               fill=TRUE,
               col.names=paste0("V", seq_len(ncol))
)
# cl
# class(cl)

## add cluster sizes for uniqued reads grouped by XhoI position (see cluster.pl)

# add new plot to current plot without cleaing the frame
par(new=TRUE)

matpoints(7:91,
          cl,
          pch=20,
          col=rgb(0,1,0, 0.5),
          xaxt="n", # don't plot x-axis
          yaxt="n" # don't plot y-axis
)
axis(4, col="green")
mtext("cluster size", side=4, line=3, col="green")

## ---- cluster size by SbfI position ----

ncol <- max(
        count.fields("all_ind_pre_SbfI_cl_size_by_pos.cl", sep=",")
)
# ncol
cl <- read.csv("all_ind_pre_SbfI_cl_size_by_pos.cl",
               sep=",",
               row.names=1,
               header=FALSE,
               fill=TRUE,
               col.names=paste0("V", seq_len(ncol))
)
# cl
# class(cl)

par(mfrow=c(1,1), bg="white")
matplot(7:39,
        cl,
        pch=20,
        col=gray(0.5, 0.5),
        xlab="position in read (1-based)",
        ylab="cluster size"
)

## ---- read pair mapping analysis ----

# source functions from John Kruschke's DBDA book
source("~/Dropbox/Kruschke_Dog_Book/DBDA2Eprograms/DBDA2E-utilities.R")
source("~/Dropbox/Kruschke_Dog_Book/DBDA2Eprograms/BernBeta.R")

# prior prob. that a read pair maps concordantly
# described by a uniform beta(theta | a=1, b=1) distribution
priorAB = c(1,1) # this is more vague than actually true, 
# but I want to let the data speak exclusively

# data from random read pairs
N = 846583 # total mapping read pairs
z = 440015 # number of read pairs mapping concordantly
posteriorAB_1 = c(z+1, (N-z)+1)

# data from SbfI site containing read pairs
N = 2184 # total mapping read pairs
z = 333 # number of read pairs mapping concordantly
posteriorAB_2 = c(z+1, (N-z)+1)

a1 = posteriorAB_1[1]
b1 = posteriorAB_1[2]
a2 = posteriorAB_2[1]
b2 = posteriorAB_2[2]

# take random sample from posterior distribution of theta for random read pairs
ThetaRandReadPairs = rbeta(10000, a1, b1)
# take random sample from posterior distribution of theta for SbfI containing read pairs
ThetaSbfIContReadPairs = rbeta(10000, a2, b2)
# now take the difference in the sampled theta values from the two posterior distributions
postDiffInTheta = abs(ThetaRandReadPairs - ThetaSbfIContReadPairs)
# plot this sample from the the posterior distribution of theta differences
plot(density(postDiffInTheta), xlim=c(0.3, 0.4),
     xlab=bquote("difference in " * theta),
     main=bquote("posterior distribution of the difference in " * theta),
     zero.line=T
)
# get info from 'density' function
densityOut = density(postDiffInTheta, n=2^10)
# turn density into prob. mass
probMass = densityOut$y / sum(densityOut$y)
# get HDI info
HDIinfo = HDIofGrid(probMassVec = probMass, credMass = 0.95) # 'HDIofGrid is function in DBDA2E-utilities.R
# get HDI limits
n = length(HDIinfo$indices)
HDI_lower = densityOut$x[HDIinfo$indices[1]]
HDI_upper = densityOut$x[HDIinfo$indices[n]]
# shade HDI area
points(densityOut$x[HDIinfo$indices], densityOut$y[HDIinfo$indices], type="h", col="skyblue")
# the following code is taken and only slightly adapted from John Kruschke's 'BernGrid.R'
# Mark vertical lines at ends of HDI (sub-)intervals:
for ( limit in c(HDI_lower, HDI_upper) ) {
        lines( c(limit, limit) , 
               c(-0.5,HDIinfo$height*sum(densityOut$y)) , # 'height' is prob. mass, turned back into density
               type="l" , 
               lty=2 , 
               lwd=1.5 
        )
        text( limit , HDIinfo$height*sum(densityOut$y) , 
              bquote(.(round(limit,3))) ,
              adj=c(0.5,-0.1) , cex=1.2 )
}
# place HDI text label
text( mean(c(HDI_upper, HDI_lower)) , HDIinfo$height*sum(densityOut$y) ,
      bquote( .(100*signif(HDIinfo$mass,3)) * "% HDI" ) ,
      #adj=c(1.1,-1.5) , 
      #cex=1.5 
)
# add label for mode of dist
modeTheta = densityOut$x[ which.max( densityOut$y ) ]
text(.32, 40, 
     labels=bquote("mode = " * .(signif(modeTheta, 3)))
     )
