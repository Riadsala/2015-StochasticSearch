library(ggplot2)
library(lme4)

d <- read.csv("FixHumanObs.txt", colClasses = c("factor", "numeric", "numeric","numeric","numeric"))
d <- subset(d, RMS==1.1)
d$difficulty <- factor(as.factor(d$beta), labels=c("rough","medium","smooth"))
names(d)[1] <- "observer"
d$ecc <- d$ecc * 1.7/100

#plt <- ggplot(d, aes(x=observer, y=nFix))  +geom_boxplot()+facet_wrap(~difficulty) 
#plt <- plt + scale_y_log10(name="number of fixations", breaks=c(1,2,4,6,8,10,20,40,60,80,100,200), minor_breaks=NULL)
#plt <- plt + theme_bw()


plt <- ggplot(d, aes(x=ecc, y=nFix, colour=observer))  +geom_smooth(method="glm", se=F, forumla=nFix~ecc*difficulty+(ecc*difficulty|observer), family="poisson")+facet_wrap(~difficulty) 
plt <- plt + scale_y_log10(name="mean number of fixations", breaks=c(1.5,2,4,6,8,10,20,40,60,80, 100), minor_breaks=NULL)
plt <- plt + scale_x_continuous(name="distance from centre of screen to target (degrees of visual angle)") + theme_bw()

ggsave("numFixObs.pdf", width=10, height = 5)

m <- glmer(data=d, nFix~ecc*beta+(ecc*beta|observer), family="poisson")

