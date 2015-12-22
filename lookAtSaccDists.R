#setwd("ScothasticSearch")
library(ggplot2)

saccDistCorner <- read.csv("saccadeDists/Corner_t=1.csv")

 pol2cart <- function(dat)
{
	names(dat) = c("r", "theta")
	dat$X = dat$r * cos(dat$theta)
	dat$Y = dat$r * sin(dat$theta)
	return(dat)
}

saccDistCorner <- pol2cart(saccDistCorner)


ggplot(saccDistCorner, aes(x=X, y=Y)) + geom_point()

