setwd("/Users/Claudius/Data_analysis/SNP-indel-calling/data")

# ---- prequel to the following chunk ----
rm(list=ls())
# depths per site per ind from 'samtools depth':
depth.table = read.table(
  "ParEry.noSEgt2.nogtQ99Cov.noDUST.3.15.noTGCAGG.ANGSD_combinedNegFisFiltered.noGtQ99GlobalCov.sorted.depths.gz", 
  header=F)
# provide column names:
names(depth.table)[1:2] = c("contig_ID", "position")
#head(depth.table)
# samtools doesn't provide column names to I assume the columns come in the order the filenames are globbed by the shell,
# i. e. in the order they were provided to samtools depth:
filenames = scan("slim.bamfile.list", what="character")
Names = gsub("Data/", "", gsub(".sorted.*", "", filenames))
names(depth.table)[3:ncol(depth.table)] = Names
#colnames(depth.table)
# save for later:
save(depth.table, file="depth.table.RData")
# get contingency table of the across individual per site depths:
global.depth = rowSums(depth.table[,3:38])
save(global.depth, file="global.depth.RData")




### ---- global-coverage-dist ----
load("global.depth.RData")
# tabulate returns a count is for each depth from 1 till max coverage:
global.depth.dist = tabulate(global.depth)
# # and the new max depth is:
# length(global.depth.dist) # which is the gloabl Q99
names(global.depth.dist) = as.character(1:length(global.depth.dist))
# plot the new global per site coverage distribution
par(mfrow=c(1,1))
barplot(global.depth.dist, xlab="across sample coverage", ylab="number of sites", 
        main="Global coverage distribution")
# average coverage per-site per-ind.:
#cov = sum(global.depth.dist*1:length(global.depth.dist))/(36*nrow(depth.table))
cov = sum(global.depth)/36/length(global.depth)



### ---- ind-cov-dist ----
load("depth.table.RData")
#names(depth.table)
par(mgp = c(4,1,0), mar = c(5,5,4,2))
z = boxplot(depth.table[,3:38], outline=F, plot=F)
bxp(z, outline=F, boxfill=c(rep("red", 18), rep("green", 18)), xaxt="n", xlab="individuals", ylab="coverage", main="Individual coverage distributions")
axis(side=1, at=1:36, labels=names(depth.table)[3:38], cex.axis=.6, las=2)




### ---- mapQ-dist ----
load("../ANGSD/Data/mapQ_ERY_afterFiltering.RData")
load("../ANGSD/Data/mapQ_ERY.RData")
load("../ANGSD/Data/mapQ_PAR_afterFiltering.RData")
load("../ANGSD/Data/mapQ_PAR.RData")
# ERY
ery.mq = as.matrix(rbind(mq.ery/sum(mq.ery), mq.after.ery/sum(mq.after.ery)))
rownames(ery.mq) = c("before filtering", "after filtering")
par(mfrow=c(1,1), mar=c(5,4,2,0))
mp = barplot(ery.mq,
             ylim=c(0,.3), 
             xlab="mapping quality score",
             ylab="proportion of reads",
             main="mapping quality distribution for ERY",
             col=c("orange", "blue"), 
             beside=TRUE,
             legend=T,
             args.legend = list(x="top")
)
axis(side=1, 
     at=mp[c(1, 10, 20, 30, 40)*2]-mp[1]/2, 
     labels=as.character(c(1, 10, 20, 30, 40))
)
# PAR
par.mq = as.matrix(rbind(mq.par/sum(mq.par), mq.after.par/sum(mq.after.par)))
rownames(par.mq) = c("before filtering", "after filtering")
par(mar=c(5,2,2,0))
mp = barplot(par.mq,
             ylim=c(0,.3), 
             xlab="mapping quality score",
             #ylab="proportion of reads",
             yaxt = "n",
             main="mapping quality distribution for PAR",
             col=c("orange", "blue"), 
             beside=TRUE,
             legend=T,
             args.legend = list(x = "top")
)
axis(side=1, 
     at=mp[c(1, 10, 20, 30, 40)*2]-mp[1]/2, 
     labels=as.character(c(1, 10, 20, 30, 40))
)
#----------------------------------



### ---- PCA ----
covar = as.matrix(read.table("EryPar.covar.noSNPcall.unknownMinor", header=F))
pca = prcomp(covar)
pc1_prop_var = (pca$sdev[1]^2)/sum((pca$sdev)^2)
pc2_prop_var = (pca$sdev[2]^2)/sum((pca$sdev)^2)
par(mfrow=c(1,1))
plot(pca$rot[,1], pca$rot[,2], xlim=c(-0.3, 0.3), 
     xlab=paste("PC1 (", signif(pc1_prop_var, digits=2)*100, "%)", sep=""), 
     ylab=paste("PC2 (", signif(pc2_prop_var, digits=2)*100, "%)", sep=""), 
     type="n",
    main="1,730,389 sites\n unknown minor allele + weighting by p(var)",
     cex.main=.9
)
points(pca$rot[,1], pca$rot[,2], 
     col=c(rep("red", 18), rep("green", 18)), 
     pch=c(rep(17, 18), rep(19, 18)), 
)
legend("topleft", 
       legend=c("ERY", "PAR"),
       pch=c(17,19),
       col=c("red", "green"),
       bty="n"
)
filenames = scan("slim.bamfile.list", what="character")
Names = gsub("Data/", "", gsub(".sorted.*", "", filenames))
# get indeces for outliers:
mu = mean(pca$rot[,1])
pc1_sd = (pca$rot[,1] - mu)^2
#Names[order(pc1_sd)]
#Names[which(pc1_sd <= sort(pc1_sd)[2])] # two closest closest inds on PC1
index = which(pc1_sd <= sort(pc1_sd)[2]) 
mu2 = mean(pca$rot[,2])
pc2_sd = (pca$rot[,2] - mu2)^2
#Names[order(pc2_sd, decreasing=TRUE)]
#Names[which(pc2_sd >= sort(pc2_sd, dec=TRUE)[3])] # three most distant inds to mean of PC2
index = c(index, which(pc2_sd >= sort(pc2_sd, dec=TRUE)[3]))
# label outliers:
text(pca$rot[,1][index], pca$rot[,2][index], labels=Names[index], cex=0.7, adj=c(0,0))



# ---- permut-global-fst-hist ----
load("perm.global.fst.RData") # creates perm.fst
load("boot.resample.Fst.by.contig.RData") # creates boot.resample.Fst.by.contig
load("bhatia.RData")
Fst.bhatia.global = sum(bhatia[["Hb.minus.Hw"]])/sum(bhatia[["Hb"]])
bhatia.CI95 = quantile(boot.resample.Fst.by.contig, probs=c(.025, .975))
perm.fst.CI95 = quantile(perm.fst, probs=c(.025, .5, .975))
#range(perm.fst)
hist(boot.resample.Fst.by.contig, 
     border="navyblue", col="navyblue", 
     freq=FALSE, 
     xlim=c(0,.4),
     main=expression(paste("global average ", F[ST])),
     xlab=expression(paste("global Bhatia's ", F[ST]))
)
hist(perm.fst, breaks=20, border="springgreen", col="springgreen", 
     add=TRUE,
     freq=FALSE
)
legend("topleft", 
       legend=c("100 random permutations of population label", 
                paste("10,000 bootstrap resamples of ", length(levels(bhatia$contig)), " contigs")),
       col=c("springgreen", "navyblue"),
       cex=.7,
       fill=c("springgreen", "navyblue")
)



# ---- fst-by-ascertainment-class ----
load("global.fst.by.ERY.ascert.RData")
load("global.fst.by.PAR.ascert.RData")
n = 18 # number of individuals in each pop
par(mar=c(5, 4, 1, 2) + 0.1)
plot(1:(n-1)/(2*n), global.fst.by.ERY.ascert,
     ylim=c(0.1,.4), ylab=expression(paste("average Bhatia's ", F[ST])),
     xlab="minor allele frequency",
     pch=17, col="red",
     type="b"
    # main=expression(paste("Allele frequency dependence of ", F[ST]))
)
lines(1:(n-1)/(2*n), global.fst.by.PAR.ascert,  
      ylim=c(0,1), ylab=expression(paste("average Bhatia's ", F[ST])),
      xlab="minor allele frequency",
      pch=19, col="green",
      type="b"
)
abline(h=Fst.bhatia.global, lwd=2, col="darkgrey")
legend("topright", 
       legend=c("ascertainment in ERY", "ascertainment in PAR", "without ascertainment"),
       col=c("red", "green", "darkgrey"),
       pch=c(17, 19, NA),
       lty=c(NA, NA, 1),
       lwd=c(NA, NA, 2),
       #fill=c("green", "red", "grey"),
       bty="n"
)
par(mar=c(5, 4, 4, 2) + 0.1)



# ---- fst-by-contig-hist ----
load("Fst.by.contig.bhatia.RData")
hist(Fst.by.contig$FST, xlab=expression(paste("average Bhatia's ", F[ST])), 
     breaks=100,
     main=bquote(paste("Distribution of ", F[ST], " from ", .(nrow(Fst.by.contig)), " contigs")),
     col=gray(.3),
     border="darkgrey"
)



# ---- unfolded-2D-SFS ----
sfs2d = scan("EryPar.unfolded.2dsfs")
sfs2d = matrix(sfs2d, nrow=37, ncol=37) # rows should be PAR, columns should be ERY
sfs2d[1,1] = 0
library(fields)
# rows in the matrix are on the x-axis, columns are on the y-axis:
image.plot(0:37, 0:37, sfs2d, xlab="PAR", ylab="ERY", main="global unfolded 2D-SFS")



# ---- folded-sfs ----
# read in folded spectra from only overlapping sites
ery.sfs = scan("../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/original/ERY/ERY.unfolded.sfs.folded")
par.sfs = scan("../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/original/PAR/PAR.unfolded.sfs.folded")
# read in 1D SFS from also non-overlapping sites
ery.sfs.nonO = scan("../ANGSD/SFS/with_ANGSD-0.917-142-ge3dbeaa/ERY/ERY.unfolded.sfs.folded")
par.sfs.nonO = scan("../ANGSD/SFS/with_ANGSD-0.917-142-ge3dbeaa/PAR/PAR.unfolded.sfs.folded")
# unfortunately, function definitions when read by knitr expire at the end
# of the code chunk
# read in functions from external file:
source("functions.R")
par(mar=c(5,4,2,0))
mp = barplot(rbind(ery.sfs[-1], ery.sfs.nonO[-1], par.sfs[-1], par.sfs.nonO[-1])/c(sum(ery.sfs[-1]), sum(ery.sfs.nonO[-1]), sum(par.sfs[-1]), sum(par.sfs.nonO[-1])), 
             names.arg=1:length(ery.sfs[-1]),
             beside=TRUE,
             col=c("red", "salmon1", "seagreen", "palegreen2"),
             xlab="minor allele count",
             ylab="proportion of segregating sites",
             main="ML folded site frequency spectrum"
)
SNM = snm(1, 18, 36)
lines(colMeans(mp), SNM/sum(SNM), pch=21, col="black", type="b")
legend("topright",
       legend=c("erythropus, only overlapping", "erythropus, including non-overlapping", "parallelus, only overlapping", "parallelus, including non-overlapping", "standard neutral model"),
       fill=c("red", "salmon1", "seagreen", "palegreen2", NA),
       border = c("black", "black", "black", "black", NA),
       pch = c(NA, NA, NA, NA, 21),
       bty = "n"
)



# ---- fit-neutral-theta ----
# get observed folded spectra from only overlapping sites:
ery.sfs = scan("../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/original/ERY/ERY.unfolded.sfs.folded")
par.sfs = scan("../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/original/PAR/PAR.unfolded.sfs.folded")
# read in functions from external file:
source("functions.R")
# unfortunately, function definitions when read by knitr expire at the end
# of the code chunk
# optimize theta
ery_thetaOpt = optimize(f, interval=c(0, sum(ery.sfs[-1])),  eta=ery.sfs[-1], n=36, maximum=FALSE, tol=0.001)
par_thetaOpt = optimize(f, interval=c(0, sum(par.sfs[-1])),  eta=par.sfs[-1], n=36, maximum=FALSE, tol=0.001)





#
#
#

# ---- folded-sfs-boot ----
# get 200 bootstrap replicates of ERY spectrum from only overlapping sites, with exhaustive search
pp = pipe("cat ../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/bootstrap/ERY/*.unfolded.sfs.folded", "r")
ery.sfs.boot = read.table(pp, header=F)
close(pp)
ery.sfs.boot = ery.sfs.boot[,-1] # discard count of monomorphic sites
# get 200 bootstrap replicates of PAR spectrum from only overlapping sites, with exhaustive search
pp = pipe("cat ../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/bootstrap/PAR/*.unfolded.sfs.folded", "r")
par.sfs.boot = read.table(pp, header=F)
close(pp)
par.sfs.boot = par.sfs.boot[,-1]

# calculate 95% bootstrap CI:
#names(ery.sfs.boot) = as.character(0:18)
ery.sfs.boot.CI95 = apply(ery.sfs.boot, 2, quantile, probs=c(0.025, 0.5, 0.975))
par.sfs.boot.CI95 = apply(par.sfs.boot, 2, quantile, probs=c(0.025, 0.5, 0.975))
#
# get optimal standard neutral model spectra
source("functions.R")
snm_ery = snm(ery_thetaOpt$min, len=18, n=36)
snm_par = snm(par_thetaOpt$min, len=18, n=36)
# need to get snm.expect now to set the proper ylim
y_max =  max(ery.sfs.boot.CI95[3,], par.sfs.boot.CI95[3,], snm_ery, snm_par)
par(mfrow=c(1,1), mar=c(5,4,2,0))
#
# ERY:
#
# plot medians of bootstraps:
plot(1:18, 
     ery.sfs.boot.CI95[2,1:18], 
     ylim=c(0, y_max),
     pch=20,
     xlab="minor allele count",
     ylab="number of sites",
     main="erythropus",
     type="b",
     col="red",
     xaxp=c(1, 18, 18-1)
)
abline(h=seq(0, y_max, 2000), lty="dashed", col="lightgrey")
arrows(1:18, ery.sfs.boot.CI95[2,1:18], 
       1:18, ery.sfs.boot.CI95[3,1:18],
       angle=90,
       length=.05,
       col="red"
)
arrows(1:18, ery.sfs.boot.CI95[2,1:18], 
       1:18, ery.sfs.boot.CI95[1,1:18],
       angle=90,
       length=.05,
       col="red"
)

#
# SNM expected:
#
points(1:18, snm_ery, pch=18, col=gray(0, 0.5), type="b")
#
legend("topright", bty="n",
       legend=c("standard neutral model"),
       pch=18,
       col="gray"
)
# # get lower and upper 95% quantiles:
# low = qpois(p=0.025, lambda=snm.expect)
# high = qpois(p=0.975, lambda=snm.expect)
# # add 95% CI bars:
# arrows(1:length(snm.expect), snm.expect, 
#        1:length(snm.expect), high,
#        angle=90,
#        length=.05
# )
# arrows(1:length(snm.expect), snm.expect, 
#        1:length(snm.expect), low,
#        angle=90,
#        length=.05
# )
#
# # plot observed spectrum:
# points(1:18, ery.sfs, pch=4, cex=1, col="red", type="b", lwd=2)
# legend(x=10, y=8200, legend="real sample", bty="n", col="red", pch=4, lwd=2, cex=.9)
# SNM expected:
# points(1:length(snm.expect), snm.expect, 
#        pch=18, col="gray", type="b"
# )
# # get lower and upper 95% quantiles:
# low = qpois(p=0.025, lambda=snm.expect)
# high = qpois(p=0.975, lambda=snm.expect)
# #
#

# PAR:
#
# plot medians of bootstraps:
par(mar=c(5, 2,2,0))
plot(1:18, 
     par.sfs.boot.CI95[2,1:18], 
     ylim=c(0, y_max),
     pch=20,
     type="b",
     xlab="minor allele count",
     ylab="",
     yaxt = "n",
     main="parallelus",
     col="seagreen",
     xaxp=c(1, 18, 18-1)
)
#axis(2, at=c(0, 5000, 10000, 15000), labels=FALSE)
abline(h=seq(0, y_max, 2000), lty="dashed", col="lightgrey")
arrows(1:18, par.sfs.boot.CI95[2,1:18], 
       1:18, par.sfs.boot.CI95[3,1:18],
       angle=90,
       length=.05,
       col="seagreen"
)
arrows(1:18, par.sfs.boot.CI95[2,1:18], 
       1:18, par.sfs.boot.CI95[1,1:18],
       angle=90,
       length=.05,
       col="seagreen"
)
#
# SNM expected:
#
points(1:18, snm_par, pch=18, col=gray(0, 0.5), type="b")
#
legend("topright", bty="n",
       legend=c("standard neutral model"),
       pch=18,
       col="gray"
)




#
#
#
#
# ---- folded-sfs-boot-exh-par ----
# PAR
par.sfs.boot.exh = read.table("PAR.FOLDED.sfs.boot.exh", header=F)
names(par.sfs.boot.exh) = as.character(0:18)
# calculate 95% CI:
par.sfs.boot.exh.CI95 = as.data.frame(t( apply(par.sfs.boot.exh, 2, quantile, probs=c(0.25, 0.5, 0.975)) ))
#
source("functions.R")
snm.expect = snm(par_thetaOpt$min, len=length(ery.sfs), n=36)
#
# plot medians of bootstraps:
plot(1:18, 
     par.sfs.boot.exh.CI95[2:19,2], 
     ylim=c(0, max(par.sfs.boot.exh.CI95[2:19,], snm.expect)),
     pch=20,
     xlab="minor allele count",
     ylab="number of sites",
     main="parallelus",
     type="b",
     col="green"
)
arrows(1:18, par.sfs.boot.exh.CI95[2:19,2], 
       1:18, par.sfs.boot.exh.CI95[2:19,3],
       angle=90,
       length=.05,
       col="green"
)
arrows(1:18, par.sfs.boot.exh.CI95[2:19,2], 
       1:18, par.sfs.boot.exh.CI95[2:19,1],
       angle=90,
       length=.05,
       col="green"
)
# 
# SNM neutral spectrum
#
# SNM expected:
points(1:length(snm.expect), snm.expect, 
       pch=18, col="gray", type="b"
)
# get lower and upper 95% quantiles:
low = qpois(p=0.025, lambda=snm.expect)
high = qpois(p=0.975, lambda=snm.expect)
# add 95% CI bars:
arrows(1:length(snm.expect), snm.expect, 
       1:length(snm.expect), high,
       angle=90,
       length=.05
)
arrows(1:length(snm.expect), snm.expect, 
       1:length(snm.expect), low,
       angle=90,
       length=.05
)
# # plot observed spectrum:
# points(1:18, par.sfs, pch=4, cex=1, col="green", type="b", lwd=2)
# legend(x=10, y=12300+500, legend="real sample", bty="n", col="green", pch=4, lwd=2, cex=.9)
legend("topright", bty="n",
       legend=c("observed", "optimal neutral fit"),
       fill=c("green", "gray"),
       border=c("green", "gray")
)










#
#
#
#
#
#

# ---- S-and-Pi ----
rm(list=ls())
# read in folded spectra from only overlapping sites
ery.sfs = scan("../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/original/ERY/ERY.unfolded.sfs.folded")
par.sfs = scan("../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/original/PAR/PAR.unfolded.sfs.folded")
# read in 1D SFS from also non-overlapping sites
ery.sfs.nonO = scan("../ANGSD/SFS/with_ANGSD-0.917-142-ge3dbeaa/ERY/ERY.unfolded.sfs.folded")
par.sfs.nonO = scan("../ANGSD/SFS/with_ANGSD-0.917-142-ge3dbeaa/PAR/PAR.unfolded.sfs.folded")
# unfortunately, function definitions when read by knitr expire at the end
# of the code chunk
# read in functions from external file:
source("functions.R")
#
# real sample SFS from overlapping sites:
#
# S:
ery.nSites = sum(ery.sfs)
S.ery = sum(ery.sfs[-1])
#--
#n = 36
#a1 = sum(sapply(1:(n-1), function(x) x^(-1)))
#theta_Watt_ery = S.ery/a1
#Nery = theta_Watt_ery/(4*3e-9*ery.nSites)
#Tery = 500000/(2*Nery) # 500ky in genetic units
#--
S_prop.ery = S.ery/ery.nSites
# pi:
pi.ery = PI( ery.sfs[-1] )
pi_sites.ery = pi.ery/ery.nSites
# S:
par.nSites = sum(par.sfs)
S.par = sum(par.sfs[-1])
#theta_Watt_par = S.par/a1
#Npar = theta_Watt_par/(4*3e-9*par.nSites)
#Tpar = 500000/(2*Npar)
#--

S_prop.par = S.par/par.nSites
pi.par = PI( par.sfs[-1] )
pi_sites.par = pi.par/par.nSites
#
# bootstrap resampled SFS from overlapping sites:
#
# get 200 bootstrap replicates of ERY spectrum from only overlapping sites, with exhaustive search
pp = pipe("cat ../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/bootstrap/ERY/*.unfolded.sfs.folded", "r")
ery.sfs.boot = read.table(pp, header=F)
close(pp)
# get 200 bootstrap replicates of PAR spectrum from only overlapping sites, with exhaustive search
pp = pipe("cat ../ANGSD/BOOTSTRAP_CONTIGS/minInd9_overlapping/SFS/bootstrap/PAR/*.unfolded.sfs.folded", "r")
par.sfs.boot = read.table(pp, header=F)
close(pp)
# number of segregating sites:
S.ery.boot = rowSums( ery.sfs.boot[,-1] )
S_prop.ery.boot = S.ery.boot/rowSums(ery.sfs.boot)
#formatC(quantile(S_prop.ery.boot, probs=c(.025, .975)), format="f", digits=4)
S.par.boot = rowSums( par.sfs.boot[,-1] )
S_prop.par.boot = S.par.boot/rowSums(par.sfs.boot)
# average number of pairwise differences:
source("functions.R")
pi.ery.boot = apply( ery.sfs.boot[,-1], 1, function(sfs) PI( sfs ) )
pi_sites.ery.boot = pi.ery.boot/rowSums(ery.sfs.boot)
#formatC(quantile(pi_sites.ery.boot, probs=c(.025, 0.5, .975)), format="f", digits=4)
pi.par.boot = apply( par.sfs.boot[,-1], 1, function(sfs) PI( sfs ) )
pi_sites.par.boot = pi.par.boot/rowSums(par.sfs.boot)
#quantile(pi_sites.par.boot, probs=c(.025, 0.5, .975))
#
# including non-overlapping sites
#
# S:
ery.nSites.nonO = sum(ery.sfs.nonO)
S.ery.nonO = sum(ery.sfs.nonO[-1])
S_prop.ery.nonO = S.ery.nonO/ery.nSites.nonO
# pi:
pi.ery.nonO = PI( ery.sfs.nonO[-1] )
pi_sites.ery.nonO = pi.ery.nonO/ery.nSites.nonO
# S:
par.nSites.nonO = sum(par.sfs.nonO)
S.par.nonO = sum(par.sfs.nonO[-1])
S_prop.par.nonO = S.par.nonO/par.nSites.nonO
pi.par.nonO = PI( par.sfs.nonO[-1] )
pi_sites.par.nonO = pi.par.nonO/par.nSites.nonO
#
# bootstrap resampled SFS from including non-overlapping sites:
#
pp = pipe("cat ../ANGSD/BOOTSTRAP_CONTIGS/including_non-overlapping/SFS/bootstrap/ERY/*.unfolded.sfs.folded", "r")
ery.sfs.boot.nonO = read.table(pp, header=FALSE)
close(pp)
#
pp = pipe("cat ../ANGSD/BOOTSTRAP_CONTIGS/including_non-overlapping/SFS/bootstrap/PAR/*.unfolded.sfs.folded", "r")
par.sfs.boot.nonO = read.table(pp, header=FALSE)
close(pp)
#
# number of segregating sites in bootstrap replicates
S.ery.boot.nonO = rowSums(ery.sfs.boot.nonO[,-1])
S_prop.ery.boot.nonO = S.ery.boot.nonO/rowSums(ery.sfs.boot.nonO)
#quantile(S_prop.ery.boot.nonO, probs=c(.025, 0.5, .975))
#
S.par.boot.nonO = rowSums(par.sfs.boot.nonO[,-1])
S_prop.par.boot.nonO = S.par.boot.nonO/rowSums(par.sfs.boot.nonO)
#quantile(S_prop.par.boot.nonO, probs=c(.025, 0.5, .975))
#
# pi from bootstrap replicates
pi.ery.boot.nonO = apply( ery.sfs.boot.nonO, 1, function(sfs) PI( sfs[-1] ) )
pi_sites.ery.boot.nonO = pi.ery.boot.nonO/rowSums(ery.sfs.boot.nonO)
#quantile(pi_sites.ery.boot.nonO, probs=c(.025, 0.5, .975))
#
pi.par.boot.nonO = apply( par.sfs.boot.nonO, 1, function(sfs) PI( sfs[-1] ) )
pi_sites.par.boot.nonO = pi.par.boot.nonO/rowSums(par.sfs.boot.nonO)
#quantile(pi_sites.par.boot.nonO, probs=c(.025, 0.5, .975))
#
# # plot S:
# plot(density(S.ery.boot/ery.nSites),
#      xlab=expression(S[prop]),
#      xlim=range(S.ery.boot/ery.nSites, S.par.boot/par.nSites),
#      main="Proportion of segregating sites")
# lines(density(S.par.boot/par.nSites), lty=2, lwd=1.5)
# points(c(S.ery/ery.nSites, S.par/par.nSites), c(0, 0), pch=3)
# legend("topright", legend=c("erythropus", "parallelus", "real sample"), 
#        lty=c(1, 2, NA), lwd=c(1, 1.5, NA), pch=c(NA, NA, 3), bty="n")
# # plot pi:
# plot(density(pi.par.boot/par.nSites),
#      xlim=range(c(pi.ery.boot/ery.nSites, pi.par.boot/par.nSites)),
#      xlab=expression(pi[site]), lty=2, lwd=1.5,
#      main="Average number of pairwise differences per nucleotide site"
# )
# lines(density(pi.ery.boot/ery.nSites), lty=1)
# points(c(pi.ery/ery.nSites, pi.par/par.nSites), c(0, 0), pch=3)
# legend("top", legend=c("erythropus", "parallelus", "real sample"), lty=c(1, 2, NA), 
#        lwd=c(1, 1.5, NA), pch=c(NA, NA, 3), bty="n")



# ---- Tajimas-D ----
# 1. calculate constant (see p. 45 in Gillespie)
# ery
n = 36
a1 = sum(sapply(1:(n-1), function(x) x^(-1)))
a2 = sum(sapply(1:(n-1), function(x) x^(-2)))
b1 = (n+1)/(3*(n-1))
b2 = 2*(n^2+n+3)/(9*n*(n-1))
c1 = b1 - (1/a1)
c2 = b2 - (n+2)/(a1*n)+a2/a1^2
C = sqrt( c1/a1*S.ery + (c2/(a1^2+a2))*S.ery*(S.ery-1) )
ery.TajimasD.global = (pi.ery - S.ery/a1)/C
# par
C = sqrt( c1/a1*S.par + (c2/(a1^2+a2))*S.par*(S.par-1) )
par.TajimasD.global = (pi.par - S.par/a1)/C
# bootstraps:
C.boot = sqrt( c1/a1*S.ery.boot + (c2/(a1^2+a2))*S.ery.boot*(S.ery.boot-1) )
ery.TajimasD.global.boot = (pi.ery.boot - S.ery.boot/a1)/C.boot
#
C.boot = sqrt( c1/a1*S.par.boot + (c2/(a1^2+a2))*S.par.boot*(S.par.boot-1) )
par.TajimasD.global.boot = (pi.par.boot - S.par.boot/a1)/C.boot
#
# including non-overlapping sites
#
# ery
C = sqrt( c1/a1*S.ery.nonO + (c2/(a1^2+a2))*S.ery.nonO*(S.ery.nonO-1) )
ery.TajimasD.global.nonO = (pi.ery.nonO - S.ery.nonO/a1)/C
# par
C = sqrt( c1/a1*S.par.nonO + (c2/(a1^2+a2))*S.par.nonO*(S.par.nonO-1) )
par.TajimasD.global.nonO = (pi.par.nonO - S.par.nonO/a1)/C
#
# bootstrap replicates:
#
# ery
C.boot.nonO = sqrt( c1/a1*S.ery.boot.nonO + (c2/(a1^2+a2))*S.ery.boot.nonO*(S.ery.boot.nonO-1) )
ery.TajimasD.global.boot.nonO = (pi.ery.boot.nonO - S.ery.boot.nonO/a1)/C.boot.nonO
#quantile(ery.TajimasD.global.boot.nonO, probs=c(0.025, 0.5, 0.975))
# par
C.boot.nonO = sqrt( c1/a1*S.par.boot.nonO + (c2/(a1^2+a2))*S.par.boot.nonO*(S.par.boot.nonO-1) )
par.TajimasD.global.boot.nonO = (pi.par.boot.nonO - S.par.boot.nonO/a1)/C.boot.nonO
#quantile(par.TajimasD.global.boot.nonO, probs=c(0.025, 0.5, 0.975))
#
# par(mar=c(5, 4, 2, 2) + 0.1)
# plot(density(ery.TajimasD.global.boot),
#      xlim=range(ery.TajimasD.global.boot, par.TajimasD.global.boot),
#      type="l",
#      xlab="Tajima's D",
#      main="global Tajima's D"
# )
# lines(density(par.TajimasD.global.boot), lty=2, lwd=1.5)
# points(c(ery.TajimasD.global, par.TajimasD.global), c(0, 0), pch=3)
# legend("top", legend=c("erythropus", "parallelus", "real sample"), lty=c(1, 2, NA), 
#        lwd=c(1, 1.5, NA), pch=c(NA, NA, 3), bty="n")







#
#
#
#
##
#
# ---- Tajimas-D-by-contig ----
load("thetas.RData")
par(mfrow=c(2,1))
hist(thetas.ery$Tajima, breaks=40, 
     main="erythropus",
     xlab="Tajima's D",
     col=gray(.3),
     border="darkgrey"
)
legend("topright", legend=paste(nrow(thetas.ery), " contigs"), bty="n")
hist(thetas.par$Tajima, breaks=40,
     main="parallelus",
     xlab="Tajima's D",
     col=gray(.3),
     border="darkgrey"
)
legend("topright", legend=paste(nrow(thetas.par), " contigs"), bty="n")
par(mfrow=c(1,1))
par(mar=c(5, 4, 4, 2) + 0.1)
