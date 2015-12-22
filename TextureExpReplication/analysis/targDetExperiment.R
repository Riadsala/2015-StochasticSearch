setwd("Documents/GIT-general-private/2014 StochasticSearch/TextureExpReplication/analysis")

library(ggplot2)
library(lme4)
library(car)
people = c(1,2,3,4,5,6,7,8,9,11)
alldat = data.frame(participantNumber=character(), trialNum=numeric(), TargetPresent=numeric(), beta=numeric(), seed=numeric(), r=numeric(), phi=numeric(), response=numeric()) 
for (person in people)
{
	dat = read.csv(paste("../results/ac_", person, "_td.txt", sep=""))
	if (person == 7)
	{
		# there's a bug in the data causing participantNumber to = 6 for particiapnt 7!
		dat$participantNumber = 7
	}
	if (person == 4)
	{
		# they got their numbers the wrong way around
		dat$respose = 1-dat$respose
	}

	names(dat)[8] = "response"
	alldat = rbind(alldat, dat)
}

visdeg = 60*tan(54.2/57)/(1920*pi)

alldat$beta <- factor(as.numeric(alldat$beta), labels =c("rough", "medium", "smooth"))
alldat$phi <- pi/180 * alldat$phi
alldat$r <- alldat$r * visdeg
alldat$x <- alldat$r*sin(alldat$phi)
alldat$y <- alldat$r*cos(alldat$phi)
alldat$x2 <- alldat$x^2
alldat$y2 <- alldat$y^2
alldat$participantNumber = as.factor(alldat$participantNumber)
alldat$TargetPresent = as.factor(alldat$TargetPresent)
alldat$rejTrial = as.numeric(alldat$response == -1)
summary((aggregate(data=alldat, rejTrial~participantNumber, FUN=mean)))
alldat = alldat[-which(alldat$response==-1),]
alldat = alldat[-which(alldat$response==2),]
accdat = aggregate(data=alldat, response~participantNumber+beta+TargetPresent, FUN=mean)
summary(aggregate(data=alldat[which(alldat$TargetPresent==0),], response~participantNumber, FUN=mean))
levels(accdat$TargetPresent) = c("target absent", "target present")
accdat$response[which(accdat$TargetPresent=="target absent")] = 1- accdat$response[which(accdat$TargetPresent=="target absent")] 


fpplt = ggplot(accdat, aes(x=participantNumber, y=response, fill=beta)) + geom_bar(stat="identity", position=position_dodge())+facet_wrap(~TargetPresent)
fpplt = fpplt + scale_y_continuous(name="accuracy", expand=c(0,0), breaks=c(seq(0,1,0.1))) + theme_bw()
fpplt = fpplt + scale_x_discrete(name="participant ID")+ scale_fill_brewer(palette="Set1")
fpplt
ggsave("targDetFalsePositive.pdf", width=10, height=4)

targPres = subset(alldat, TargetPresent==1 & response!=-1 & response!=2)

targPres2 <- aggregate(response ~ participantNumber+ beta+r, targPres, FUN="mean")
tmp <- aggregate(response ~ beta+r, targPres, FUN="var")
targPres2$stderr = sqrt(tmp$response/9)
plt1 <- ggplot(targPres2, aes(x=r, y=response, colour=beta, shape=beta)) + geom_smooth(method="glm", family="binomial", se=F) + geom_point(position = position_jitter(w = 0.15, h = 0.0))
plt1 <- plt1 + scale_y_continuous(name="p(target detected)")  +scale_x_continuous(name="eccentricity (degrees of visual angle)") + theme_bw()+ scale_colour_brewer(palette="Set1")
plt1
ggsave("targDetModelEccOnly.pdf", width=10, height=5)

plt2 <- ggplot(targPres2, aes(x=r, y=response, colour=beta)) + geom_smooth(method="glm", family="binomial", se=F) + geom_point(position = position_jitter(w = 0.15, h = 0.0))
plt2 <- plt2 + scale_y_continuous(name="p(target detected)")  +scale_x_continuous(name="eccentricity (degrees of visual angle)") + theme_bw()
plt2 = plt2 + facet_wrap(~participantNumber)
ggsave("eccAccPltByPerson.pdf")

m2 <- glm(data=targPres, response ~ beta * (x2+y2), family="binomial")

targPres3 <- aggregate(response ~ participantNumber+ beta+r+phi, targPres, FUN="mean")
targPres3$phi = targPres3$phi * 180/phi
library(dplyr)
targPres3 = filter(targPres3, phi %in% c(0, 180))

# contour map
gd <- as.data.frame(expand.grid(x=seq(-10,10,0.1),y=seq(-10,10,0.1)))
n <- length(gd$x)
cD <- data.frame(x=rep(gd$x,3), y=rep(gd$y,3), beta=factor(rep(c(1,2,3), c(n,n,n)),labels=c("rough","medium","smooth")))
cD$x2 <- cD$x^2
cD$y2 <- cD$y^2
cD$z = predict(m2, cD, type="response")

cplt <- ggplot(cD, aes(x=x, y=y, z=z)) + coord_fixed() +  stat_contour(binwidth=0.1,size=0.75,aes(colour = ..level..)) 
cplt = cplt + facet_wrap(~beta) + theme_bw()
ggsave("aggtargdet2.pdf", width=10, height=5)


mm = glmer(data=targPres, response ~ beta * x2+y2 + (beta|participantNumber), family="binomial")
