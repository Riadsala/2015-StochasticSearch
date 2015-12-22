library(ggplot2)
library(mvtnorm)
# Start the clock!
ptm <- proc.time()

maxFix = 1

# parameters 
dprimeE = 20 

dprimeI_atFix = 5/dmvnorm(c(0,0),sigma=array(c(2,0,0,2), dim=c(2,2)))

# set up stimuli 
# i will count item locations
n = 3
item = data.frame(ID= -0.5, x = rep(1:n, n), y=rep(1:n,each=n))
item$ID[2] = 0.5
nItems 			= nrow(item)
prior 	 		= rep(1/nItems, nItems)


possFixLocations.x = item$x
possFixLocations.y = item$y

dprimeI <- function(fix)
{
	# calculates dprime for each item location from every fixation location
	out = array(dim=c(nrow(fix), nrow(item)))
	for (t in 1:nrow(fix))
	{
		out[t,] = dprimeI_atFix * dmvnorm(t(rbind(item$x, item$y)), mean=c(fix$x[t], fix$y[t]), sigma=array(c(2,0,0,2), dim=c(2,2)))
	}
	return(out)
}

g <- function(k,t,T)
{
	if (T>1)
	{
		topfrac = dprimeE^2 * dpI2[t,]
    	botfrac = dprimeE^2 + colSums(dpI2[1:T,])
	}
	else
	{
		topfrac = dprimeE^2 * dpI2[t,]
   		botfrac = dprimeE^2 + (dpI2[1:T,])
	}
    return(topfrac/botfrac)

}

calc_pT <- function()
{
	# calculate prob for each item location
	print(T)
	sum_g_W_over_T = 0
	for (t in 1:T)
	{
		sum_g_W_over_T = sum_g_W_over_T + g(k,t,T)*W[t,]
	}
	pT = prior * exp(sum_g_W_over_T)	
	# normalise
	pT = pT / sum(pT)
	calc_pT <- pT
}



sigma <- function(T)
{
	# equation A50 - used in computation of where to look next
	# computes answers for all target locations and returns vector
	# note, eq A50 uses T+1, but as sigma function is only ever called 
	# as sigma(T+1), the +1 is taken care of before we get to this function. 
	topfrac = dprimeE^2 + colSums(dpI2)
	botfrac =  dpI2[T,] * (topfrac-dpI2[T,])
	return(sqrt(topfrac/botfrac))
}

probCorrect <- function(i,k,T)
{
	# equation A39 - used in computation of where to look next
		# used in eq A39
	print(k)
	
	h <- function(i,j,T)
	{
		# equation A36 - used in computation of where to look next
		h = log(prior[i]/prior[j])
		
		for (t in 1:(T-1))
		{
			gT = g(k=ktest,t,T)
			h = h + gT[i]*W[t,i] - gT[j]*W[t,j]	
		}
		#print(paste('h =', h))
		return(h)
	}

	u <- function(i,T)
	{
		# equation A48 and A49 - used in computation of where to look next
		targconstant = rep(-0.5, nrow(item))
		targconstant[i] = 0.5

		if (T==1)
		{
			topfrac = (dpI2[1:(T),] * (W[1:(T),]-targconstant))
			botfrac = dprimeE^2 + (dpI2[1:(T),])
		}
		else
		{
			topfrac = colSums(dpI2[1:(T),] * (W[1:(T),]-targconstant))
			botfrac = dprimeE^2 + colSums(dpI2[1:(T),])
		}
		return((topfrac/botfrac)+targconstant)
	}

	phi <- function(x)
	{
		return(1/sqrt(2*pi) * exp(-0.5*x^2))
	}

	Phi <- function(x) 
	{
		output = array()
		# this is confusing... x is a vector due to numerical integration sfuff I think
		for (z in 1:length(x))
		{
			num.int = integrate(phi, lower = -10, upper = x[z])
			output[z] = num.int$value
		}
		return(output)
	}

	integrand <- function(inputvar)
	{
		a = 1
		
		for (j in 1:nItems)
		{
			if (i!=j)
			{
				topfrac = g0[i] * (sig[i]*inputvar + uj[j]) + h(i,j,T+1)
				botfrac = g0[j] * sig[j]
				x = (topfrac/botfrac - uj[j]/sig[j])
				#print(paste('topfrac',  topfrac))
				a = a * Phi(x)
				rm(x, topfrac, botfrac)
			}
		}
		
		return((1/sig[i]) * phi(inputvar)*a)
	}
	uj = u(i,T)
	num.int = integrate(integrand, lower=-10, upper=10)
	#print(paste('propCorrect', num.int$value))
	return(num.int$value)
}



# add external noise
X = rnorm(length(item$ID), mean=0, sd=1/dprimeE)
# inital fixation at origin
k = data.frame(x=0, y=0)
# W is the "template response" at each position and fixation
W = array(dim=c(maxFix+1, length(item$ID)))
# probabilty of target pT at each position and fixation
pT = array(dim=c(maxFix, nItems))


# for each fixation
for (T in 1:maxFix )
{
	print(paste("Fixation", T))
	# generate internal noise
	dpI = dprimeI(k)
	dpI2 = dpI^2
	N = rnorm(nItems, mean=0, sd=1/dpI[T,])
	# calculate template responses
	W[T,] = item$ID + N + X
	rm(N)
	# calculate probability at each location i
	pT[T,] = calc_pT()

    print('working out where to look next')
	# now where do we look next?
	probCorrectForPosFix = data.frame(x=numeric(), y=numeric, p=numeric)
	for (x in possFixLocations.x)
	{
		for (y in possFixLocations.y)
		{
			#print(paste(x,y))
			ktest = rbind(k, data.frame(x=x, y=y))
			dpI = dprimeI(ktest)
			dpI2 = dpI^2
			# pre compute g, sigma, etc to save a bit of time
			g0 = g(k=ktest, T+1, T+1)
			sig = sigma(T+1)
			p = 0
			for (i in 1:nItems)
			{
				p = p + pT[T, i] * probCorrect(i, ktest, T)
			}
			probCorrectForPosFix = rbind(probCorrectForPosFix, data.frame(x=x,y=y,p=p))
		}
	}

	argmax = probCorrectForPosFix[which(probCorrectForPosFix$p==max(probCorrectForPosFix$p)),]
	print(argmax)
	k = rbind(k, data.frame(x=argmax$x,y=argmax$y))

	item$pT = pT[T,]


#p2 = ggplot(probCorrectForPosFix, aes(x=x, y=y, z=p))  + geom_contour(aes(colour=..level..), bins=30) + scale_gradient_brewer(type="seq", palette=2)
 p2 = ggplot(probCorrectForPosFix, aes(x=x, y=y, z=p))  + geom_tile(aes(fill = p)) + stat_contour()+   scale_fill_gradient(low = "brown", high = "white")
p2 = p2 + geom_point(data=item, aes(x=item$x, item$y, size=pT, colour=ID)) + scale_size(range = c(1, 10))
p2
p2 = p2 + geom_path(data=k, aes(x=x, y=y), colour="orange") + geom_point(data=k, aes(x=x, y=y), colour="orange")
ggsave(paste("plt", T, ".pdf", sep=""))
}

# Stop the clock
proc.time() - ptm
 # plot


