# check levels of IOR and compare simulations to humans
library(ggplot2)

source("getSaccStats.R")
source("getSaccByRegion.R")


humDat <- read.csv("../clarke2009data/clarke2009Fixations.txt", sep=" ")
humDat$strategy <- "human"
names(humDat)[2] <- "trial"


humDat$x[humDat$x<1] = NaN
humDat$y[humDat$y<1] = NaN
humDat$x[humDat$x>1024] = NaN
humDat$y[humDat$y>1024] = NaN

## do human data
ctr = 0;
d2 <- array()
d1 <- array()
theta1 <- array()
theta2 <- array()
qdat = data.frame(qfrom=character(), x=numeric(), y=numeric())
for  (o in 1:7)
{
	obsDat = subset(humDat, Obs==o)

	for (t in levels(obsDat$trial)) { 
		trialDat = subset(obsDat, trial==t)	
		saccStats = getSaccStats(trialDat)
		qs = getSaccByRegion(trialDat)
     	qdat = rbind(qdat, qs)

		theta1 = c(theta1, saccStats$theta1)
		theta2 = c(theta2, saccStats$theta2)
		d1 = c(d1, saccStats$d1)
		d2 = c(d2, saccStats$d2)
		rm(trialDat, saccStats)
	}
	rm(obsDat)
}



