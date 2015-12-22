
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
people = c(4,5)
source("basicAnalysis.R")


fixDat = data.frame(person=numeric(), trial=numeric(), tp=numeric(), resp=numeric(), beta=character(), n=numeric(), x=numeric(), y=numeric())
saccDat = data.frame(person=numeric(), trial=numeric(), amp=numeric(), ang=numeric(), fromX=numeric(), fromY=numeric())

for (person in people)
{
	print(person)
	asc = readLines(paste("results/ac_", person, "_vs.asc", sep=""))
	asc = strsplit(asc, "\t")

	trialStarts = grep("SYNCTIME", asc)
	trialEnds = grep("TRIAL_OVER", asc)

	nTrials = length(trialStarts)

	as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

	n0 = 21

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
	 		trialDat$n = 1:length(trialDat$x)
	 		trialDat$x[which(trialDat$x<1 | trialDat$x>1024)] = NaN
	 		trialDat$y[which(trialDat$y<1 | trialDat$y>1024)] = NaN
	 		saccDat = rbind(saccDat, saccStats(n, trialDat))
	 		fixDat = rbind(fixDat, trialDat)
	 		rm(trialDat)
	 		# get distance from target to closest fixation
	 		#targX = 
	 	}
}
}
fixDat$beta = as.factor(fixDat$beta)

# only look at TP trials for which we have a response
tpTrials = fixDat[which(fixDat$tp==1 & fixDat$resp==1),]

numFixDat= aggregate(data=tpTrials, x~person+trial+beta, FUN=length)
numFixDat$person = as.factor(numFixDat$person)
names(numFixDat)[4] = "numFix"

write.csv(numFixDat, "numberFixations.txt", quote=F, row.names=F)



library(ggplot2)

plt = ggplot(numFixDat, aes(x=beta, y=numFix, fill=person)) + geom_boxplot()
plt = plt + scale_y_continuous(name="number of fixations", trans = "log2", breaks=c(2,4,8,16,32,64,128))
plt
ggsave("circularPerf.pdf", height=5, width=10)

###
# let us try and make a hotspot map
fixDat = fixDat[-which(fixDat$n==1),]
fixDatQ  = fixDat
fixDatQ$x = round(fixDatQ$x/50)
fixDatQ$y = round(fixDatQ$y/50)
dat = aggregate(data=fixDatQ, person ~ x + y, FUN=length)
names(dat)[3] = "z"

plt2 = ggplot(dat, aes(x=x,y=y,z=(z))) + stat_contour(geom="polygon", aes(fill=..level..))
ggsave("hotspot.pdf")

fixDat = fixDat[complete.cases(fixDat),]
write.csv(data.frame(person=fixDat$person, x=fixDat$x,y=fixDat$y), "allFixations.txt", quote=F, row.names=F)
