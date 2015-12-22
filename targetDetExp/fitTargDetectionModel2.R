setwd("Documents/GIT-general-private/2014 StochasticSearch/targetDetExp/")
library(ggplot2)
library(lme4)
library(car)

dat <- read.csv("targDetResponses.txt", colClasses=c("factor", "factor", "numeric", "numeric", "numeric"))
dat$beta <- factor(as.numeric(dat$beta), labels =c("rough", "medium", "smooth"))
dat$obs <- factor(dat$obs, labels = c("observer 1", "observer 2"))
dat$r <- dat$r * 0.84 / 50
dat$phi <- pi/180 * dat$phi
dat$x <- dat$r*sin(dat$phi)
dat$y <- dat$r*cos(dat$phi)
dat$x2 <- dat$x^2
dat$y2 <- dat$y^2


# first look at target absent trials

targAb <- subset(dat, r==Inf)
aggregate(100*resp~obs+beta, targAb, FUN="mean")

targPres <- subset(dat, r!=Inf)




targPres2 <- aggregate(resp~obs+beta+r, targPres, FUN="mean")
plt1 <- ggplot(targPres2, aes(x=r, y=resp, colour=beta)) + geom_smooth(method="glm", family="binomial", se=F) + geom_point() + facet_wrap(~obs)
plt1 <- plt1 + scale_y_continuous(name="p(target detected)")  +scale_x_continuous(name="eccentricity (degrees of visual angle)") + theme_bw()
ggsave("targDetModelEccOnly.pdf", width=10, height=5)


m <- glm(data=targPres, resp ~ beta * x2 + y2, family="binomial")
summary(m)

m <- glm(data=targPres, resp ~ beta * (x2 + y2), family="binomial")
summary(m)

# contour map
gd <- as.data.frame(expand.grid(x=seq(-10,10,0.1),y=seq(-10,10,0.1)))
n <- length(gd$x)
cD <- data.frame(x=rep(gd$x,3), y=rep(gd$y,3), beta=factor(rep(c(1,2,3), c(n,n,n)),labels=c("rough","medium","smooth")))
cD$r = with(cD, sqrt(x^2+y^2))
cD$x2 <- cD$x^2
cD$y2 <- cD$y^2
cD$z = predict(m, cD, type="response")

cplt <- ggplot(cD, aes(x, y, z=z)) + coord_fixed() +  stat_contour(binwidth=0.1,size=0.75,aes(colour = ..level..)) + facet_wrap(~beta) 
cplt <- cplt + theme_bw()+scale_colour_gradient(name="probability of detection", limits=c(0, 1), low="red", high="green") + xlim(c(-8.6,8.6)) + ylim(c(-8.6,8.6))
ggsave("newTDcontourplot.pdf", width=10, height=5)