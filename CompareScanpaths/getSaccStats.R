

getSaccStats <- function(trialDat)
{
	d1 <- array()
	d2 <- array()
	theta1 <- array()
	theta2 <- array()

	n = max(trialDat$f)	
	if (n>1)
	{
		ctr <- 0
		for (f in 2:n)	 	
		{
		   ctr = ctr + 1
		   # compute distance from fixation f to fixation f-1
		   dx <- trialDat$x[f]-trialDat$x[f-1]
		   dy <- trialDat$y[f]-trialDat$y[f-1] 
		   d1[ctr] <- (sqrt(dx^2 + dy^2))
		   theta1[ctr] <- 180/pi * atan2(dy,dx)	
		}
		# check to see if we have enough data to compute two-step stats for this trial	
		if (n>2) 	
		{
			 ctr <- 0
			 for (f in 3:n)	 	
			 {
			   	ctr = ctr + 1
			   	# compute distance from fixation f to fixation f-2
			   	ddx = trialDat$x[f]-trialDat$x[f-2]
			   	ddy = trialDat$y[f]-trialDat$y[f-2]
			   	d2[ctr] = (sqrt(ddx^2+ddy^2))
			   	
			   	# calculate relative angle
			   	dx1 <- trialDat$x[f]-trialDat$x[f-1]
		   		dy1 <- trialDat$y[f]-trialDat$y[f-1]
		   		phi1 <- 180/pi * atan2(dy1,dx1)	 
				dx2 <- trialDat$x[f-1]-trialDat$x[f-2]
		   		dy2 <- trialDat$y[f-1]-trialDat$y[f-2] 
		   		phi2 <- 180/pi * atan2(dy2,dx2)
		   		theta2[ctr] <- phi1-phi2

				
			}
		}
	}
	out = list(theta1=theta1, theta2=theta2, d1=d1, d2=d2)
}