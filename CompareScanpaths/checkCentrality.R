# check how central fixations are
# compare datasets with LLH


library(mvtnorm)
library(ggplot2)
simDat <- read.csv("../matlabSearchSim/scanpaths.txt")
simDat = subset(simDat, n>1)
humDat <- read.csv("../saccadeDists/clarke2009Fixations.txt", sep=" ")
humDat$strategy <- "human"
fixDat <- data.frame(strategy=c(simDat$strategy, humDat$strategy), n=c(simDat$n, humDat$f), x=c(simDat$x, humDat$x), y=c(simDat$y, humDat$y))
levels(fixDat$strategy) = c("optimal", "stochastic", "human")
fixDat$x <- (fixDat$x-512)/512
fixDat$y <- (fixDat$y-512)/512
fixDat$x[fixDat$x < (-1)] = NaN
fixDat$x[fixDat$x > (1)] = NaN
fixDat$y[fixDat$y < (-1)] = NaN
fixDat$y[fixDat$y > (1)] = NaN 

# reshape to x and y are like factors
n <- length(fixDat$x)
dat = data.frame(strategy = rep(fixDat$strategy,2), d = c(fixDat$x, fixDat$y), axis = rep(c("x","y"),c(n,n)))



ggplot(dat, aes(x=d, fill=strategy)) + geom_density(alpha=0.3) + facet_wrap(~axis)


sig2 = 0.22
nu = 0.45
Sig <- c(sig2, 0, 0, nu*sig2)
dim(Sig) <- c(2,2)

ss <- subset(fixDat, strategy=="optimal")
p = dmvnorm(as.matrix(cbind(ss$x, ss$y)), c(0,0), Sig, log=T)
# randomly sample N points to assess how well dist fits data
pr <- sample(na.omit(p), 1000)
logLikOptimal = sum(pr)

ss <- subset(fixDat, strategy=="random")
p = dmvnorm(as.matrix(cbind(ss$x, ss$y)), c(0,0), Sig, log=T)
# randomly sample N points to assess how well dist fits data
pr <- sample(na.omit(p), 1000)
logLikRandom = sum(pr)

#


