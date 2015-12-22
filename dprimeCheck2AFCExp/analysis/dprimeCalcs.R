library(psyphy)
library(dplyr)
library(ggplot2)
dat2afc = read.csv('parsedResults.txt', header=T, sep=' ')
dat2afc$beta = as.factor(dat2afc$beta)
levels(dat2afc$beta) = c("rough", "medium", "smooth")
dat2afc$p[which(dat2afc$p==1)] = 0.99
dat2afc$phi = round(dat2afc$phi * 180/pi)
dat2afc$ecc = 100*round(dat2afc$ecc/100)
dat2afc$phi[which(dat2afc$phi==360)] = 0
dat2afc = filter(dat2afc, ecc!=0)
dat2afc = aggregate(p ~ beta+ ecc,dat2afc, FUN=mean)

for (i in 1:nrow(dat2afc))
{
	dat2afc$dprime[i] = dprime.mAFC(dat2afc$p[i],2)
}

visdeg = 60*tan(54.2/57)/(1920*pi)

dat2afc$ecc = visdeg * dat2afc$ecc
plt = ggplot(dat2afc, aes(x=ecc, y=dprime, colour=beta))
plt = plt + geom_point() + geom_smooth(method=lm, se=F)
plt = plt + theme_bw() + theme(legend.position="none")

plt = plt + scale_y_continuous(name="d'", limit=c(0,4))
plt= plt + scale_x_continuous(name="target eccentricity (visual degrees)")
ggsave("dprime.pdf", width=4, height=4)
# plt

## we want to load yes/no data.

people = c(1:9,11)
alldat = data.frame(participantNumber=character(), trialNum=numeric(), TargetPresent=numeric(), beta=numeric(), seed=numeric(), r=numeric(), phi=numeric(), response=numeric()) 
for (person in people)
{
	dat = read.csv(paste("../../TextureExpReplication/results/ac_", person, "_td.txt", sep=""))
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

# visdeg = 60*tan(54.2/57)/(1920*pi)

# alldat$beta <-  factor(as.numeric(alldat$beta), labels =c("rough", "medium", "smooth"))
# alldat$phi <- pi/180 * alldat$phi
names(alldat)[6]="ecc"
alldat$participantNumber = as.factor(alldat$participantNumber)
alldat$TargetPresent = as.factor(alldat$TargetPresent)
alldat$beta = as.factor(alldat$beta)
# alldat$rejTrial = as.numeric(alldat$response == -1)
# summary((aggregate(data=alldat, rejTrial~participantNumber, FUN=mean)))

alldat = filter(alldat,  TargetPresent==1, response!=2, participantNumber==11)#
# alldat = filter(alldat, phi%in%c(0,180))
datYesNo = aggregate(data=alldat, response~beta+ecc, FUN=mean)


fp= aggregate(data=dat, response~TargetPresent+beta, FUN=mean)
datYesNo$response[which(datYesNo$response == 1)]=0.999
datYesNo$dprime[datYesNo$beta==1.6] = qnorm(datYesNo$response[datYesNo$beta==1.6])-qnorm(fp$response[fp$TargetPresent==0 & fp$beta==1.6])
datYesNo$dprime[datYesNo$beta==1.65] = qnorm(datYesNo$response[datYesNo$beta==1.65])-qnorm(fp$response[fp$TargetPresent==0 & fp$beta==1.65])
datYesNo$dprime[datYesNo$beta==1.7] = qnorm(datYesNo$response[datYesNo$beta==1.7])-qnorm(fp$response[fp$TargetPresent==0 & fp$beta==1.7])
levels(datYesNo$beta) = c("rough", "medium", "smooth")


datYesNo$ecc = visdeg * datYesNo$ecc
plt = ggplot(datYesNo, aes(x=ecc, y=dprime, colour=beta))
plt = plt + geom_point() + geom_smooth(method=lm, se=F)
plt = plt + theme_bw()
plt = plt + scale_y_continuous(name="d'", limit=c(0,4))
plt= plt + scale_x_continuous(name="target eccentricity (visual degrees)")
ggsave("YesNo_AllPhi.pdf", width=5, height=4)



