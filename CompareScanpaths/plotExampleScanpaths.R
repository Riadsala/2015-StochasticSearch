
library(ggplot2)

source("getSaccStats.R")

simDat <- read.csv("../searchsim/exScanpathsForPlot.txt")

levels(simDat$strategy) = c("optimal", "stochastic")
humDat <- read.csv("../clarke2009data/clarke2009Fixations.txt", sep=" ")

# find a trial which had 20 fixations
idx = which(humDat$f==20)
trialDat = humDat[(idx[1]-21):idx[1],]

scanpaths = data.frame(observer = factor(c(simDat$sim, 3*trialDat$Obs), labels=c("optimal", "stochastic", "human")), f = c(simDat$t, trialDat$f),x = c(simDat$x, trialDat$x),y = c(simDat$y, trialDat$y) )

scanpaths$observer = factor(scanpaths$observer, levels(scanpaths$observer)[c(3,2,1)])

splt <- ggplot(scanpaths, aes(x=x, y=y, colour=observer))+ geom_path() + geom_point() + facet_wrap(~observer) + theme_bw() + coord_fixed()
splt <- splt + scale_x_continuous(limits = c(1, 1024), breaks=NA) + scale_y_continuous(limits = c(1, 1024), breaks = NA)
ggsave("exScanpaths.pdf", width=10, height=5)