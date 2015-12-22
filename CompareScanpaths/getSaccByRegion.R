

getSaccByRegion <- function(trialDat)
{
	
	N = 1024
	Q = 3

    q = N/Q


	n = max(trialDat$f)	

qfrom = array()
x= array()
y = array()


	if (n>1)
	{
		ctr <- 0
		for (f in 2:n)	 	
		{
		   ctr = ctr + 1
		   qx <- ceiling(trialDat$x[f-1]/q)
		   qy <- ceiling(trialDat$x[f-1])
		   qfrom[ctr] <- c(as.character(qx), as.character(qy))
		   x[ctr] <- trialDat$x[f]
		   y[ctr] <- trialDat$y[f]

		}
		


	}
	saccades <- rbind(saccades, data.frame(qfrom=qfrom, x=x, y=y))
	out = saccades
}