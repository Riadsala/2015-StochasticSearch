library(ggplot2)


saccStats <- function(t, fix)
{
 	nfix = nrow(fix)
 	saccs = data.frame(person=numeric(), trial=numeric(), amp=numeric(), dir=numeric(), fromX=numeric(), fromY=numeric())
 	if (nfix>1)
 	{ 	
	 	for (f in 1:(nfix-1))
	 	{
	 		d = sqrt((fix$x[f+1]-fix$x[f])^2 + (fix$y[f+1]-fix$y[f])^2) 
	 		phi = atan2(fix$x[f+1]-fix$x[f], fix$y[f+1]-fix$y[f])
	 		saccs = rbind(saccs, data.frame(person=person, trial=t, amp=d, dir=phi, fromX=fix$x[f], fromY=fix$y[f],toX=fix$x[f+1], toY=fix$y[f+1]))
	 	}
 	}
 	return(saccs)
}

# pixels to visual degrees facotr
toVD = 1# 0.84 / 50

# parse eye-tracking data
people = c(1,2,3,4,5,6,7,9)
source("basicAnalysis.R")

trialInfo$correct = trialInfo$RT>0
trialInfo$correct = as.numeric(trialInfo$correct == trialInfo$TargetPresent)
aggregate(data=trialInfo, correct ~ TargetPresent, FUN=mean)

targPres = trialInfo[which(trialInfo$TargetPresent==1),]
targPres = targPres[-which(targPres$RT==-1),]
targPres$RT = log(targPres$RT)
targPresAg = aggregate(data=targPres, RT~beta, FUN=mean)
targPresAg2 = aggregate(data=targPres, RT~beta, FUN=sd)
targPresAg$sder = targPresAg2$RT/sqrt(length(people))
rm(targPresAg2)
pltRT = ggplot(targPresAg, aes(x=beta, y=RT)) + geom_line()
pltRT = pltRT + geom_errorbar(aes(x=beta,y=RT, ymin=RT-1.96*sder, ymax=RT+1.96*sder))
pltRT = pltRT + scale_y_continuous(name="mean log RT") + theme_bw()
ggsave("meanLogRT.pdf")


library(lme4)
library(car)
m.rt = lmer(data=targPres, RT ~ beta + (beta|participantNumber))
Anova(m.rt)
fixDat = data.frame(person=numeric(), trial=numeric(), tp=numeric(), resp=numeric(), beta=character(), n=numeric(), x=numeric(), y=numeric())
saccDat = data.frame(person=numeric(), trial=numeric(), amp=numeric(), ang=numeric(), fromX=numeric(), fromY=numeric())

for (person in people)
{
	print(person)
	if (person == 8)
	{
asc = readLines(paste("../results/subj8/ac_8_1.asc", sep=""))
	}
	else
	{asc = readLines(paste("../results/ac_", person, "_vs.asc", sep=""))}
	asc = strsplit(asc, "\t")

	trialStarts = grep("SYNCTIME", asc)
	trialEnds = grep("TRIAL_OVER", asc)

	nTrials = length(trialStarts)

	as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

	if (sum(person==c(1,2))==1) {n0 = 1}
	else {n0=21}

	for (n in n0:nTrials)
	{
	 	trial = asc[trialStarts[n]:trialEnds[n]]
	 	fixationLines = grep("EFIX", trial)
	 	if (length(fixationLines)>0)
	 	{
		 	fixations = as.data.frame(matrix(unlist(trial[fixationLines]), byrow=T, ncol=6))
	 		
	 		beta = trialInfo$beta[which(trialInfo$participantNumber==person & trialInfo$trialNum==(n-n0+1))]
	 		tp = trialInfo$TargetPresent[which(trialInfo$participantNumber==person & trialInfo$trialNum==(n-n0+1))]
	 		resp = (trialInfo$RT[which(trialInfo$participantNumber==person & trialInfo$trialNum==(n-n0+1))]>0)
	 		trialDat = data.frame(person=rep(person,length(resp)), trial=n, tp=tp, resp=resp, beta=beta, 
	 			x=as.numeric.factor(fixations$V4), y=as.numeric.factor(fixations$V5))
	 		
	 		# convert to stimulus coordinates
	 		trialDat$x = trialDat$x - (1920-1024)/2
	 		trialDat$y = trialDat$y - (1080-1024)/2
	 		trialDat$x[which(trialDat$x<1 | trialDat$x>1024)] = NaN
	 		trialDat$y[which(trialDat$y<1 | trialDat$y>1024)] = NaN
	 		trialDat$n = 1:length(trialDat$x)
	 		saccDat = rbind(saccDat, saccStats(n, trialDat))
	 		fixDat = rbind(fixDat, trialDat)
	 		rm(trialDat)
	 		# get distance from target to closest fixation
	 		#targX = 
	 	}
	}
}
fixDat$beta = as.factor(fixDat$beta)

saveRDS(fixDat, "humanFixations.txt")

# only look at TP trials for which we have a response
tpTrials = fixDat[which(fixDat$tp==1 & fixDat$resp==1),]

numFixDat= aggregate(data=tpTrials, x~person+trial+beta, FUN=length)
numFixDat$person = as.factor(numFixDat$person)
names(numFixDat)[4] = "numFix"

write.csv(numFixDat, "numberFixations.txt", quote=F, row.names=F)

targPresAg2 = aggregate(data=targPres, RT~beta, FUN=sd)
targPresAg$sder = targPresAg2$RT/sqrt(length(people))

pltRT = ggplot(targPresAg, aes(x=beta, y=RT)) + geom_line()
pltRT = pltRT + geom_errorbar(aes(ymin=RT-1.96*sder, ymax=RT+1.96*sder))
pltRT = pltRT + geom_errorbar(aes(x=beta,y=RT, ymin=RT-1.96*sder, ymax=RT+1.96*sder))
pltRT = pltRT + scale_y_continuous(name="mean log RT") + theme_bw()
ggsave("meanLogRT.pdf", height=5, width=5)


plt = ggplot(numFixDat, aes(x=beta, y=numFix, fill=person)) + geom_boxplot()
plt = plt + scale_y_continuous(name="number of fixations", trans = "log2", breaks=c(2,4,8,16,32,64,128))
plt
ggsave("newVSexp.pdf", height=5, width=10)

taTrials = trialInfo[which(trialInfo$TargetPresent==0),]
taTrials$resp = taTrials$RT == -1
tagg = aggregate(data=taTrials, resp~participantNumber, FUN=mean)


theta = hist(180*saccDat$dir/pi,45)
saccAngle = data.frame(theta=theta$mids%%360, count=theta$density)
saccAmp = data.frame(amp=saccDat$amp/60)
write.table(saccDat[complete.cases(saccDat),], "saccDat.txt", quote=F, row.names=F, col.names=F)

ggplot(saccAmp, aes(x=amp))+ geom_histogram() + scale_x_continuous(name="saccadic amplitude", limits=c(0,20))+ scale_y_continuous(name="number of saccades")+theme_bw()
ggsave("saccAmp.pdf", width=6, height=4)

roseplt <- ggplot(saccAngle, aes(x=theta, y=count)) + geom_bar(width=10, stat="identity") 
roseplt <- roseplt + scale_x_continuous(name=" ", limits = c(0, 360),breaks = c(0, 45, 90, 135, 180, 225, 270, 315))+scale_y_continuous(name=" ",breaks=NULL)+ coord_polar(start=0, direction=-1)+theme_bw()
ggsave("roseplot.pdf", width=6, height=4)


 fixDat = fixDat[-which(fixDat$n==1),]
fixDat = fixDat[complete.cases(fixDat),]
write.csv(data.frame(person=fixDat$person, x=fixDat$x,y=fixDat$y), "allFixations.txt", quote=F, row.names=F)

