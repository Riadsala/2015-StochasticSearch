# check levels of IOR and compare simulations to humans
library(ggplot2)

source("getSaccStats.R")


humDat <- read.csv("../clarke2009data/clarke2009Fixations.txt", sep=" ")
humDat$strategy <- "human"
humDat$Obs <- as.factor(humDat$Obs)
names(humDat)[2] = "trial"
humDat$trial = as.factor(humDat$trial)


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
		

		theta1 = c(theta1, saccStats$theta1)
		theta2 = c(theta2, saccStats$theta2)
		d1 = c(d1, saccStats$d1)
		d2 = c(d2, saccStats$d2)
		rm(trialDat, saccStats)
	}
	rm(obsDat)
}

theta = hist(theta1,45)
saccAngle = data.frame(theta=theta$mids%%360, count=theta$density)
saccAmp = data.frame(amp=d1/60)

ggplot(saccAmp, aes(x=amp))+ geom_histogram() + scale_x_continuous(name="saccadic amplitude", limits=c(0,20))+ scale_y_continuous(name="number of saccades")+theme_bw()
ggsave("saccAmp.pdf", width=5, height=5)



roseplt <- ggplot(saccAngle, aes(x=theta, y=count)) + geom_bar(width=10, stat="identity") 
roseplt <- roseplt + scale_x_continuous(name=" ", limits = c(0, 360),breaks = c(0, 45, 90, 135, 180, 225, 270, 315))+scale_y_continuous(name=" ",breaks=NULL)+ coord_polar(start=3*pi/2, direction=-1)+theme_bw()
ggsave("roseplot.pdf", width=5, height=5)


