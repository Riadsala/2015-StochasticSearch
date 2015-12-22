library(ggplot2)

maxFix = 20

# parameters 
dprimeE = 2 #sqrt(c^2/(alpha/en))
dprimeI_atFix = 2

# set up stimuli 
# i will count item locations
itemIDs 		= c(0.5, -0.5, -0.5)
nItems 			= length(itemIDs)
itemLocations 	= array(dim=c(2,nItems), c(-1,0, 1,0, 5,0))
prior 	 		= rep(1/nItems, nItems)

dprimeI <- function(t, i=NULL, fix=NULL)
{
	# will use fixation t from fixations list, unless we explicitly pass a fixation
	if (is.null(fix))
	{
		fix = fixations[[t]]
	}
	# for now, lets assume a linear fall-off with distance
	if (is.null(i))
	{
		# haven't specifed an item, so calculate dprimeI for all items
		distFromFix = sqrt((itemLocations[1,]-fix[1])^2+(itemLocations[2,]-fix[2])^2)
	}
	else
	{
		# only calcualte for i-th item
		distFromFix = sqrt((itemLocations[1,i]-fix[1])^2+(itemLocations[2,i]-fix[2])^2)	
	}
	out = dprimeI_atFix/distFromFix
	return(out)
}

g <- function(t)
{
	topfrac = dprimeE^2 * dprimeI(t)^2
    botfrac = dprimeE^2 + sumdprimeI2_1toT(t)
    return(topfrac/botfrac)
}

sumdprimeI2_1toT <- function(F)
{
	sumdprimeI2_1toT = 0
	for (t in 1:F)
	{
    	sumdprimeI2_1toT = sumdprimeI2_1toT + dprimeI(t)^2
    }
    return(sumdprimeI2_1toT)
}

calc_pT <- function()
{
	# calculate prob for each item location
	sum_g_W_over_T = 0
	for (t in 1:T)
	{
		sum_g_W_over_T = sum_g_W_over_T + g(t)*W[t,]
	}
	pT[T,] = prior * exp(sum_g_W_over_T)

	
	# normalise
	pT[T,] = pT[T,] / sum(pT[T,])
	calc_pT <- pT
}

sigma <- function(potentialFix)
{
	topfrac = dprimeE^2
	for (t in 1:T)
    {
    	topfrac = topfrac + dprimeI(t,i)^2
    }
    botfrac = dprime(fix=potentialFix, i=i)^2 * (topfrac)
    topfrac = topfrac + dprime(potentialFix, itemLocations[,i])^2
    sigma = topfrac/botfrac
}


possFixLocatios = list(x = -5:5, y=-5:5)

# add external noise
X = rnorm(length(itemIDs), mean=0, sd=1/dprimeE)
# inital fixation at origin
fixations = list(c(0,0))

# W is the "template response" at each position and fixation
W = array(dim=c(maxFix, length(itemIDs)))
# probabilty of target pT at each position and fixation
pT = array(dim=c(maxFix, length(itemIDs)))

# for each fixation
for (T in 1:maxFix )
{
	print(T)
	# generate internal noise
	dPrimeIforCurrentFix = dprimeI(T)	
	N = rnorm(length(itemIDs), mean=0, sd=1/dPrimeIforCurrentFix)
	# calculate template responses
	W[T,] = itemIDs + N + X
	rm(N)
	# calculate probability at each location i
	pT = calc_pT()

	# now where do we look next?
	for (x in possFixLocatios$x)
	{
		for (y in possFixLocatios$y)
		{
			
		}
	}
	fixations = c(fixations, list(c(0,0)))

}
 # plot
dat1 = data.frame(item = rep(c(1,2,3), c(maxFix, maxFix, maxFix)), f = rep(1:maxFix, nItems), pT = c(pT[,1],pT[,2],pT[,3]))
p1 <- ggplot(dat1, aes(x=f, y=pT, colour=as.factor(item))) + geom_path()
dat2 = data.frame(item = rep(c(1,2,3), c(maxFix, maxFix, maxFix)), f = rep(1:maxFix, nItems), W = c(W[,1],W[,2],W[,3]))
p2 <- ggplot(dat2, aes(x=f, y=W, colour=as.factor(item))) + geom_path()