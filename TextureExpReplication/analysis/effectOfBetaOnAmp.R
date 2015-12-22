setwd("Documents/GIT-general-private/2014 StochasticSearch/TextureExpReplication/analysis")

# check levels of IOR and compare simulations to humans
library(ggplot2)
library(dplyr)

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

humDat = readRDS("humanFixations.txt")
humDat$person <- as.factor(humDat$person)
humDat$trial = as.factor(humDat$trial)

## do human data
ctr = 0;
saccAmps = data.frame(person=numeric(), beta=numeric(), amp=numeric())
theta1 <- array()
theta2 <- array()
for  (o in 1:7)
{
	for (b in levels(humDat$beta))
	{
		obsDat = filter(humDat, person==o,  beta==b)

		for (t in levels(obsDat$trial)) { 
			trialDat = subset(obsDat, trial==t)	
			saccStats = getSaccStats(trialDat);
			theta1 = c(theta1, saccStats$theta1)
			theta2 = c(theta2, saccStats$theta2)
			n = length(saccStats$d1)
			saccAmps = rbind(saccAmps, data.frame(person=rep(o,n), beta=rep(b,n), amp=saccStats$d1))
			rm(trialDat, saccStats)
		}
	}
	rm(obsDat)
}

saccAmps$amp = saccAmps$amp/60 

pltAmp <- ggplot(saccAmps, aes(x=amp, fill=beta))+ geom_density(alpha=0.3) 
pltAmp <- pltAmp + scale_x_continuous(name="displacement (degrees of visual angle)")  + theme_bw()
ggsave("onetwosaccamp.pdf", width=10, height=4)
levels(angDat$stat) = c("direction", "relative direction")
pltAng <- ggplot(angDat, aes(x=x, fill=strategy))+ geom_density(alpha=0.3)+ facet_wrap(~stat, scale="free")
pltAng <- pltAng + scale_x_continuous(name="saccade angle (degrees)", breaks=c(90, 180, 270))+ theme_bw()
ggsave("saccAngles.pdf", width=10, height=4)








