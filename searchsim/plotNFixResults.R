setwd("Documents/GIT-general-private/2014 StochasticSearch/searchsim/")

library(ggplot2)
simdat <- read.table("searchSimResults.txt", header = TRUE,  sep = ",",  stringsAsFactors = FALSE)
p1 <- ggplot(simdat, aes(y=log(nfix), x=beta, colour=strat)) + geom_point() + geom_smooth(method=lm) + facet_wrap(~ecc)

model <- glm(data=simdat, nfix ~ beta*ecc*strat, family="poisson")


obsdat <- read.table("FixHumanObs.txt", header = TRUE,  sep = ",",  stringsAsFactors = FALSE)
obsdat$Obs <- as.factor(obsdat$Obs)



# combine into one datset
names(obsdat)[1] = "Observer"
names(simdat)[1] = "Observer"
names(simdat)[5] = "nFix"


# subset multiple factors at once
# NOTE: subtracting 1 from simulation fixations as human observers do not need to fixate target 
# to detect it (while hte model is hard coded to fixate target if it spots it)
dat <- data.frame(
	strategy = c(rep("human", length(obsdat$Observer)), simdat$Observer),
	Observer = c(as.factor(obsdat$Observer), simdat$Observer), 
	beta=c(obsdat$beta, simdat$beta), 
	ecc=c(obsdat$ecc, simdat$ecc), 
	nFix=c(obsdat$nFix, simdat$nFix-1))

library(dplyr)
dat$beta <- factor(as.factor(dat$beta), labels=c("rough", "medium", "smooth"))
dat$ecc <- factor(as.factor(dat$ecc/60), labels=c("1.67", "3.75", "5.83"))
simdat$strategy = factor(simdat$strategy, levels(simdat$strategy)[c(1,3,2)])
# dat = dat[which(dat$strategy!="optimal"),]
p2 <- ggplot(filter(dat, strategy!='optimal'), aes(x=ecc, y=nFix, fill=strategy, cond=Observer), width=500, height=600) 
p2 <- p2 + geom_boxplot() + facet_wrap(~beta, nrow=3) + scale_y_continuous(name="number of fixations", trans = "log2")+scale_x_discrete(name="target eccentricity (degrees of visual angle)")
p2 <- p2 + theme_bw() + scale_fill_brewer(type="seq", palette=3)
p2
ggsave("numFixHumanModelclarke2009.pdf", width=6, height=6)


agDat = aggregate(data=dat, nFix ~ Observer+beta+ecc+strategy, FUN=median)

p1 <- ggplot(agDat, aes(x=ecc, y=nFix, group=Observer, colour=strategy)) 
p1 <- p1 + geom_smooth(method="lm", se=F) + facet_wrap(~beta, nrow=1) + scale_y_continuous(breaks=c(1,3,5,7,9,11))+scale_x_discrete(name="target eccentricity (degrees of visual angle)")
p1 <- p1 + theme_bw()
p1
ggsave("numFixHumanModelclarke2009lines.pdf", width=8, height=4)
agDat


n = length(simdat$nFix)
stochasticFix = simdat$nFix[seq(1,n,2)]
optimalFix = simdat$nFix[seq(2,n,2)]
t.test(stochasticFix, optimalFix, paired=T, alternative="greater")
library(car)
library(lme4)
dat$ecc = as.numeric(dat$ecc)

modelDat = dat[-which(dat$strategy=="human"),]
 mm = glm(data=modelDat, nFix ~ beta * ecc * strategy, family="poisson",)