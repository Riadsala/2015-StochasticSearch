setwd("Documents/GIT-general-private/2014 StochasticSearch/TextureExpReplication/modelling")
library(dplyr)
library(ggplot2)
simdat <- read.table("searchSimResults.txt", header = TRUE,  sep = ",",  stringsAsFactors = FALSE)
p1 <- ggplot(simdat, aes(y=log(nfix), x=beta, colour=strat)) + geom_point() + geom_smooth(method=lm) + facet_wrap(~ecc)
names(simdat)[1] = "strategy"
simdat$strategy = as.factor(simdat$strategy)
#model <- glm(data=simdat, nfix ~ beta*strat, family="poisson")


obsdat <- read.table("../analysis/numberFixations.txt", header = TRUE,  sep = ",",  stringsAsFactors = FALSE)




# combine into one datset
names(obsdat)[1] = "Observer"
names(simdat)[1] = "Observer"
names(simdat)[3] = "nFix"


# subset multiple factors at once
# NOTE: subtracting 1 from simulation fixations as human o bservers do not need to fixate target 
# to detect it (while hte model is hard coded to fixate target if it spots it)
dat <- data.frame(
	strategy = c(rep("human", length(obsdat$Observer)), simdat$Observer),
	Observer = c(as.factor(obsdat$Observer), simdat$Observer), 
	beta=c(obsdat$beta, simdat$beta), 
	nFix=c(obsdat$numFix, simdat$nFix))

dat$beta <- factor(as.factor(dat$beta), labels=c("rough", "medium", "smooth"))
dat$Observer = as.factor(dat$Observer)

levels(dat$strategy) = c("optimal", "stochastic", "human")
dat$strategy = factor(dat$strategy,levels(dat$strategy)[c(3,2,1)])
levels(dat$Observer) = seq(1,11)
dat$Observer[which(dat$strategy=="stochastic")] = 10
dat$Observer[which(dat$strategy=="optimal")] = 11



p2 <- ggplot(filter(dat, strategy!='optimal'), aes(x=beta, y=nFix, fill=strategy, cond=Observer)) 
p2 <- p2 + geom_boxplot() + scale_y_continuous(name="number of fixations" ,trans = "log2", breaks=c(1,2,4,8,16,32,64,128))+scale_x_discrete(name="target eccentricity (degrees of visual angle)")
p2 <- p2 + theme_bw()  + scale_fill_brewer(type="seq", palette=3)
p2
ggsave("numFixHumanModel.pdf", width=8, height=5)
summary(dat)


humanDat = filter(dat, strategy=='human')

stochasticFix = simdat$nFix[simdat$Observer=="stochastic"]
optimalFix = simdat$nFix[simdat$Observer=="optimal"]
t.test(stochasticFix, optimalFix, paired=T, alternative="greater")

 m = glmer(data=filter(dat, strategy!='optimal'), nFix ~ beta * strategy + (1+beta|Observer), family='poisson')
 summary(m)