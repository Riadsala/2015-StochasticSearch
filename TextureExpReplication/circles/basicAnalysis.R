## Check non-eyetracking results


trialInfo = data.frame(participantNumber=numeric(), trialNum=numeric(), TargetPresent=numeric(), betea=numeric(), seed=numeric(), x=numeric(), y=numeric(), RT=numeric())
for (person in people)
{
	dat = read.csv(paste("results/ac_", person, "_vs.txt", sep=""))
	# if (person>2)
	# {
	# 	dat = dat[-which(dat$trialNum<21),]
	# }
	trialInfo = rbind(trialInfo, dat)

	catchCorrect = mean(dat$RT[which(dat$TargetPresent==0)]==-1)
	print(paste("proportion of good catch trials:", catchCorrect))

	targetFound =  mean(dat$RT[which(dat$TargetPresent==1)]>0)
	print(paste("proportion of target's found: ", targetFound))
}


write.csv(format(trialInfo, digits=3), "trialInfo.txt", quote=F, row.names=F)
