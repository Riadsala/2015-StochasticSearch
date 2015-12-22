# check levels of IOR and compare simulations to humans
library(ggplot2)
setwd("Documents/GIT-general-private/2014 StochasticSearch/TextureExpReplication/analysis")

getSaccStats <- function(trialDat)
{
	d1 <- array()
	d2 <- array()
	theta1 <- array()
	theta2 <- array()

	n = max(trialDat$n)	
	if (n>1)
	{
		ctr <- 0
		for (f in 2:n)	 	
		{
		   ctr = ctr + 1
		   # compute distance from fixation f to fixation f-1
		   dx <- trialDat$x[f]-trialDat$x[f-1]
		   dy <- trialDat$y[f]-trialDat$y[f-1] 
		   d1[ctr] <- (sqrt(dx^2 + dy^2))
		   theta1[ctr] <- 180/pi * atan2(dy,dx)	
		}
		# check to see if we have enough data to compute two-step stats for this trial	
		if (n>2) 	
		{
			 ctr <- 0
			 for (f in 3:n)	 	
			 {
			   	ctr = ctr + 1
			   	# compute distance from fixation f to fixation f-2
			   	ddx = trialDat$x[f]-trialDat$x[f-2]
			   	ddy = trialDat$y[f]-trialDat$y[f-2]
			   	d2[ctr] = (sqrt(ddx^2+ddy^2))
			   	
			   	# calculate relative angle
			   	dx1 <- trialDat$x[f]-trialDat$x[f-1]
		   		dy1 <- trialDat$y[f]-trialDat$y[f-1]
		   		phi1 <- 180/pi * atan2(dy1,dx1)	 
				dx2 <- trialDat$x[f-1]-trialDat$x[f-2]
		   		dy2 <- trialDat$y[f-1]-trialDat$y[f-2] 
		   		phi2 <- 180/pi * atan2(dy2,dx2)
		   		theta2[ctr] <- phi1-phi2

				
			}
		}
	}
	out = list(theta1=theta1, theta2=theta2, d1=d1, d2=d2)
}

simDat <- read.csv("../modelling/scanpaths.txt")
names(simDat) = c("strategy", "trial", "beta", "n", "x", "y")
simDat$trial <- as.numeric(simDat$trial)
levels(simDat$strategy) = c("optimal", "stochastic")
humDat = readRDS("humanFixations.txt")
humDat$strategy <- "human"
humDat$person <- as.factor(humDat$person)
fixDat <- data.frame(strategy=c(simDat$strat, humDat$strategy), obs = c(rep(NaN,length(simDat$n)), humDat$person), trial = c(simDat$trial, humDat$trial), n=c(simDat$n, humDat$n), x=c(simDat$x, humDat$x), y=c(simDat$y, humDat$y))
levels(fixDat$strategy) = c("optimal", "stochastic", "human")
fixDat$trial <- as.factor(fixDat$trial)
rm(simDat, humDat)


## do stochastic model 
modelDat = subset(fixDat, strategy=="stochastic")
d2 <- array()
d1 <- array()
theta1 <- array()
theta2 <- array()
for (t in levels(modelDat$trial)){
	trialDat = subset(modelDat, trial==t)	
	saccStats = getSaccStats(trialDat);
	theta1 = c(theta1, saccStats$theta1)
	theta2 = c(theta2, saccStats$theta2)
	d1 = c(d1, saccStats$d1)
	d2 = c(d2, saccStats$d2)
	rm(trialDat, saccStats)
}
saccAmpS = data.frame(strategy="stochastic", 
	stat=c(rep("one saccade",length(d1)), rep("two saccades",length(d2))), x=c(d1, d2))

saccAngS = data.frame(strategy="stochastic", 
	stat=c(rep("theta1", length(theta1)), rep("theta2", length(theta2))), x=c(theta1, theta2))


rm(d1, d2, theta1, theta2)

# ## do optimal model
# modelDat = subset(fixDat, strategy=="optimal")
# ctr = 0;
# d2 <- array()
# d1 <- array()
# theta1 <- array()
# theta2 <- array()
# for (t in levels(modelDat$trial)){
# 	trialDat = subset(modelDat, trial==t)	
# 	saccStats = getSaccStats(trialDat);
# 	theta1 = c(theta1, saccStats$theta1)
# 	theta2 = c(theta2, saccStats$theta2)
# 	d1 = c(d1, saccStats$d1)
# 	d2 = c(d2, saccStats$d2)
# 	rm(trialDat, saccStats)
# }
# saccAmpO = data.frame(strategy="optimal",
# 	stat=c(rep("one saccade",length(d1)), rep("two saccades",length(d2))), x=c(d1, d2))

# saccAngO = data.frame(strategy="optimal", 
# 	stat=c(rep("theta1", length(theta1)), rep("theta2", length(theta2))), x=c(theta1, theta2))

# rm(d1, d2, theta1, theta2)

## do human data
ctr = 0;
d2 <- array()
d1 <- array()
theta1 <- array()
theta2 <- array()
for  (o in 1:7)
{
	obsDat = subset(fixDat, strategy=="human" & obs==o)

	for (t in levels(obsDat$trial)) { 
		trialDat = subset(obsDat, trial==t)	
		saccStats = getSaccStats(trialDat);
		theta1 = c(theta1, saccStats$theta1)
		theta2 = c(theta2, saccStats$theta2)
		d1 = c(d1, saccStats$d1)
		d2 = c(d2, saccStats$d2)
		rm(trialDat, saccStats)
	}
	rm(obsDat)
}
saccAmpH = data.frame(strategy="human", 
	stat=c(rep("one saccade",length(d1)), rep("two saccades",length(d2))), x=c(d1, d2))

saccAngH = data.frame(strategy="human", 
	stat=c(rep("theta1", length(theta1)), rep("theta2", length(theta2))), x=c(theta1, theta2))

rm(d1, d2, theta1, theta2)


# ampDat = rbind(saccAmpO, saccAmpS, saccAmpH)
# angDat = rbind(saccAngO, saccAngS, saccAngH)

ampDat = rbind(saccAmpS, saccAmpH)
angDat = rbind(saccAngS, saccAngH)

ampDat$x[which(ampDat$x>sqrt(2*512^2))] = NaN
ampDat$x = ampDat$x / 60
angDat$x <- angDat$x %% 360
ampDat$strategy = factor(ampDat$strategy, levels(ampDat$strategy)[c(3,2,1)])
angDat$strategy = factor(angDat$strategy, levels(angDat$strategy)[c(3,2,1)])


library(dplyr)
pltAmp <- ggplot(ampDat, aes(x=x, fill=strategy))+ geom_density(alpha=0.3) + facet_grid(strategy~stat, scale="free")
pltAmp <- pltAmp + scale_x_continuous(name="displacement (degrees of visual angle)")  + theme_bw()
ggsave("onetwosaccamp.pdf", width=10, height=4)
levels(angDat$stat) = c("direction", "relative direction")
pltAng <- ggplot(angDat, aes(x=x, fill=strategy))+ geom_density(alpha=0.3)+ facet_grid(strategy~stat, scale="free")
pltAng <- pltAng + scale_x_continuous(name="saccade angle (degrees)", breaks=c(90, 180, 270))+ theme_bw()
ggsave("saccAngles.pdf", width=10, height=4)








