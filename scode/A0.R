# Calibration and tmrca estimation project script A.0
# Created Saturday, September 8 2012 - Sebastian Duchene
# Modifications September 10 2012
# Need to implement models 

# Modification September 24 2012
#		This modification implements different possible rate distributions to simulate the data
#			one can choose the options "lognormal" or "exponential". If the clock is exponential, the standard deviation can be set to 0 as it will not be used in the fucntion
##################
##################
##################
##################
# Script A.0 should simulate ultrametric trees under a yule (pure-birth) speciation process. Trees are scaled for a specific root height, and then branches are multiplied by an evolutionary rate as random lognormal variable, to produce a substitution-in-branches tree.
# The user should be able to modify the root height for the trees (as a function or fixed value), and the number of tips (also as a function or a fixed value). Parameters for the rate, mean and sd can also be specified.
# Note that the the lambda (speciation rate), is always the number of tips/tree length. So this value cannot be toggled with if the tips are fixed.

##############
# sim.dataset function to simulate dataset. Returns a list, element [[1]] is the tree, element [[2]] is the chornogram,and element [[3]] is the alignment.

require(phangorn)
require(geiger)
require(laser)
require(diversitree)
## PLEASE NOTE THAT THE MEAN AND SD ARE IN THE LOG SPACE
sim.dataset <- function(seq.length, n.tax, mean.rate, sd.rate, rate.dist, root.height){

	sim.dat <- list()
	
	tr <- sim.bdtree(b = 1, d = 0, stop= "taxa",  n = n.tax)
	
	
	n.dep <- node.depth(tr)
	n.dep <- n.dep[n.tax+1:length(n.dep)]
	n.dep <- n.dep[!(is.na(n.dep))]
		
	scaling.fact <- root.height/max(branching.times(tr))

	lens <- abs(tr$edge.length)
	tr$edge.length <- lens
	lens <- lens*scaling.fact
	tr$edge.length <- lens

	tr.2 <- tr


#mean.rate <- log(rate)
#sd.rate <- log(sd)	
# The tree has been scaled to the given root height. Now the tree will be scaled using a random lognormal variable for the rate.
if(rate.dist=="lognormal"){
	branch.lens <- vector()
	for(i in 1:length(lens)){
		branch.lens[i] <- lens[i]*rlnorm(1, mean.rate, sd.rate) #rlnorm yields numbers in real space, no no need to e(1)^
		}
		
	tr.2$edge.length <- branch.lens
	}else if(rate.dist=="exponential"){
			branch.lens <- vector()
			for(j in 1:length(lens)){
				branch.lens[j] <- lens[j]*rexp(1, 1/mean.rate)
				}
	tr.2$edge.length <- branch.lens
}
# The branch lengths are now in substitutions per site according to the given rate of evolution

	#######################
	# We can now simulate the sequences
# ! Pending models and other simulation parameters. Check Paradis's book.	
	secs.tr <- as.DNAbin(simSeq(tr.2, l=seq.length))
	
	
	sim.dat[[1]] <- tr.2
	sim.dat[[2]] <- tr
	sim.dat[[3]] <- secs.tr
	sim.dat[[4]] <- data.frame(rate=mean.rate, sdrate=sd.rate, rHeight=root.height)
	
	names(sim.dat) <- c("tree_subt", "tree_time", "sequences", "sim_params")

return(sim.dat)
}	
	

