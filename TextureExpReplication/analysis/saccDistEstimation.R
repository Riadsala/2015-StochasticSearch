source("parseEyeTrack.R")

Q = 32

# quantise (x,y)'s 

saccDat[,5:8] = round(saccDat[,5:8]) %% Q + 1

A = aggregate(data=saccDat, trial ~ as.factor(fromX) + as.factor(fromY) + as.factor(toX) + as.factor(toY), FUN=length)
