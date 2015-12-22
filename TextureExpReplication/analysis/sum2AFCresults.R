
library(psyphy)
library(ggplot2)
people = c(1,2,3,4,5,6,7,8,9)
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
alldat$r <- alldat$r #* visdeg
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


sim2afc = aggregate(data=alldat, response~beta+x2+y2+TargetPresent, FUN=mean)

sim2afc$acc = sim2afc$response + (1-sim2afc$response)*0.5
sim2afc$acc[which(sim2afc$acc==1)] = 0.99
for (i in 1:nrow(sim2afc))
{
sim2afc$dprime[i] = dprime.mAFC(sim2afc$acc[i],2)

}

plt = ggplot(sim2afc, aes(x=, y=dprime)) + geom_point() + facet_wrap(~beta)
plt = plt + geom_smooth(method=lm, colour="red") + theme_bw()
ggsave("dprime_est.pdf")

summary(lm(data=sim2afc, dprime~x2*beta+y2))