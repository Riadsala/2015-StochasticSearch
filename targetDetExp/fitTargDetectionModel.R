# setwd("Desktop/StochasticSearch")
library(reshape2)
library(ggplot2)

dat1 <- read.table("acc1.csv", sep = ",")
dat2 <- read.table("acc2.csv", sep = ",")
dat1$ecc = 50*1:9
dat2$ecc = 50*1:9
dat1 <- melt(dat1, id.vars="ecc")
dat2 <- melt(dat2, id.vars="ecc")
dat = rbind(dat1, dat2)
dat$obs <- as.factor(c(rep("A", 27), rep("B", 27)))
names(dat) <- c("ecc", "beta", "prob", "obs")
levels(dat$beta) <- c("1.6", "1.65", "1.7")
dat$ecc <- dat$ecc / 60 # approx vis degrees. .. fix properly later
# fitGaussian function to data
targetDetectionModel <- nls(prob ~ a[beta] * exp(-(ecc/z[beta])^k[beta]), data=dat, start = list(a=c(.5,.5,.5), z=c(4,4,4),k=c(2,2,2)), lower=c(z=0), upper= c(a=c(1,1,1), z=c(10,10,10),k=c(5,5,5)), algorithm="port", nls.control(maxiter = 5000))
p <- ggplot(dat, aes(x=ecc, y=prob, colour=beta, shape=obs)) + geom_point() + geom_line(aes(y=fitted(targetDetectionModel)))


# now look at angle also!

dat <- read.csv("targetDetExp/rawAccDat.txt")
dat$obs <- as.factor(dat$obs)
dat$beta <- as.factor(dat$beta)
dat$phi = 45 * dat$phi - 45
dat$phi <- dat$phi * pi / 180
dat$r <- 5/6 * dat$r

dat$X = dat$r * cos(dat$phi)
dat$Y = dat$r * sin(dat$phi)


targetDetectionModel <- nls(acc ~ a[beta] * exp(-(X/zx[beta])^2-(Y/zy[beta])^2), data=dat, start = list(a=c(.5,.5,.5), zx=c(4,4,4), zy=c(4,4,4)), lower=c(z=0), upper= c(a=c(1,1,1), zx=c(10,10,10), zy=c(10,10,10)), algorithm="port", nls.control(maxiter = 5000))


dat$modelFit = predict(targetDetectionModel, dat)

p <- ggplot(dat, aes(x=phi, y=acc, colour=obs))+ geom_point() + facet_grid(beta~r)
p <- p  + geom_line(aes(y=modelFit))