# define function to optimise:
f = function(theta, eta, n){ 
  #
  # theta: parameter to optimize
  # eta: will get the observed SFS
  # n: number of gene copies sampled per diploid locus
  #
  # returns sum of squared deviations (residuals) between observed counts
  # and a candidate model
  #
  sumofsq = numeric(length(ery.sfs))
  i = seq(1, length(eta))
  delta = ifelse(i == n-i, 1, 0)
  expected = theta * (1/i + 1/(n-i)) / (1 + delta)
  sumofsq = sum( (expected - eta)^2 )
  return(sumofsq)
}
#
# define a function that returns expected counts given a theta under the assumption of
# the standard neutral model:
snm = function(theta, len=0, n=36){
  #
  # theta: should be optimized theta
  # len: should be the length of the folded SFS, i. e. highest freq. class
  # n: should be number of gene copies sampled per locus, i. e. 2*N for diploid loci
  #
  # returns expected SFS
  #
  i = seq(1, len)
  delta = ifelse(i == n-i, 1, 0)
  expected = theta * (1/i + 1/(n-i)) / (1 + delta)
  return(expected)
}

PI = function(sfs){
  n.half = 18
  n = 36
  1/(n*(n-1)/2) * sum( sapply(1:n.half, function(i) i*(n-i)*sfs[i]) )
}

#### ---- fit-neutral-theta ----
#x=10
# get observed folded spectra:
#ery.sfs = scan("ERY.FOLDED.sfs")
#par.sfs = scan("PAR.FOLDED.sfs")
#ery.sfs = ery.sfs[-1]
#par.sfs = par.sfs[-1]
# read in functions from external file:
#source("functions.R")
# unfortunately, function definitions cannot be included directly into
# this R script, as it leads to compilation errors when compiled with Knitr
# and Latex
# optimize theta
#ery_thetaOpt = optimize(f, interval=c(0, sum(ery.sfs)),  eta=ery.sfs, n=36, maximum=FALSE, tol=0.001)
#par_thetaOpt = optimize(f, interval=c(0, sum(par.sfs)),  eta=par.sfs, n=36, maximum=FALSE, tol=0.001)