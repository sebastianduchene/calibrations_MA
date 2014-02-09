# Calibration and tmrca estimation project script A.1

# MODIFICATIONS IN PROGRESS FOR FIXED NUMBER OF MULTIPLE CALIBRATIONS

# Created Saturday, September 8 2012 - Sebastian Duchene
# Modified so that the shallowest node to estimate is not the shallowest node, but that with 10% the root height

# Notes: This script should have two parts. The first should take the sim.dataset output and select a node at random. It then reports the tmrca of the node in the chornogram and the descending taxa -This is the calibration. Then it reports the tmrca for the shallowest node (!=0), its descending taxa, and those of the median node.

# In an other function, the output from the first should be used to create an xml file. The xml file will have the calibration, an explicit estimate of the shallowest and median nodes, and a fixed tree. The logfile should include the posterior, date estimates for root, calibrated node

select.calibration <- function(sim.dat, mult.cal, n.cal ,exclude.root){
	#List to report the results
	calib.sim <- list()

	#Saving general variables first
	br.times <- branching.times(sim.dat[[2]])
	br.times.m <- as.matrix(br.times)
	n.nodes <- sim.dat[[2]]$Nnode
	n.tax <- length(sim.dat[[2]]$tip.label)
	root.node <- as.numeric(rownames(br.times.m)[br.times.m==max(br.times.m)])


	#Select median point in the tree.
	#Note that this is a node closest to the MIDPOINT in the TREE. Not the median node depth
	#Note that the calibration can be the midpoint of the tree, this is reported on the output, and can be changed

	dist.nodes.to.mid.tree <- abs(br.times.m-(max(br.times.m)/2))
	med.node <- as.numeric(rownames(br.times.m)[(which(dist.nodes.to.mid.tree==min(dist.nodes.to.mid.tree)))][1]) #'which' can return several values, choose only one
	med.time <- br.times.m[rownames(br.times.m)==med.node]
	med.tax <- tips(sim.dat[[2]], med.node)

		#This can help if one actually wants the median node value
			#which(abs(br.times[-cal.node]-median(br.times[-cal.node]))==min(abs(br.times[-cal.node]-median(br.times[-cal.node]))), arr.ind=T)

	#Select the shallowest node in the tree that is different to 0
#	min.node <- as.numeric(rownames(br.times.m)[br.times==min(br.times.m[br.times.m!=0])])
#	min.time <- br.times.m[rownames(br.times.m)==min.node]
#	min.tax <- tips(sim.dat[[2]], min.node)
	dist.to.10 <- abs(br.times.m-(max(br.times.m)*0.1))
	min.node <- as.numeric(rownames(br.times.m)[ (which(dist.to.10==min(dist.to.10)))][1])
	min.time <- br.times.m[rownames(br.times.m)==min.node]
	min.tax <- tips(sim.dat[[2]], min.node)





	interest.tax <- list()
	interest.tax[[1]] <- min.tax
	interest.tax[[2]] <- med.tax

	#Select a calibration. The while loop makes sure the calibration is not set to 0, as could happen with very shallow nodes
        if(class(mult.cal) == "logical"){
	#Modifications for multiple calibrations
	#First we choose a number of calibrations = number of nodes-
            if(mult.cal==T && exclude.root==F){

		num.cal<- sample(1:(n.tax-3), 1)
		cal.nodes <- sample(as.numeric(rownames(br.times.m))[-(c(which(((n.tax+1):(2*n.tax-1))==min.node), which(((n.tax+1):(2*n.tax-1))==med.node)))] , num.cal )

		}else if(mult.cal==F && exclude.root==F){

			num.cal <- 1
			cal.nodes <- sample(as.numeric(rownames(br.times.m))[-(c(which(((n.tax+1):(2*n.tax-1))==min.node), which(((n.tax+1):(2*n.tax-1))==med.node)))] , num.cal )

		}else if(mult.cal==T && exclude.root==T){

			num.cal<- sample(1:(n.tax-4), 1)
			cal.nodes <- sample(as.numeric(rownames(br.times.m))[-(c( which( ((n.tax+1):(2*n.tax-1))==root.node),which(((n.tax+1):(2*n.tax-1))==min.node), which(((n.tax+1):(2*n.tax-1))==med.node)))] , num.cal )

		}else if(mult.cal==F && exclude.root==T){

			num.cal <- 1
			cal.nodes<- sample(as.numeric(rownames(br.times.m))[-(c( which( ((n.tax+1):(2*n.tax-1))==root.node),which(((n.tax+1):(2*n.tax-1))==min.node), which(((n.tax+1):(2*n.tax-1))==med.node)))] , num.cal )

		}
    }else if(class(mult.cal) == "character"){

# ADDITIONAL CODE TO SET A FIXED NUMBER OF CALIBRATIONS
        if(mult.cal == "fixed"){
            num.cal = n.cal
            cal.nodes<- sample(as.numeric(rownames(br.times.m))[-(c( which( ((n.tax+1):(2*n.tax-1))==root.node),which(((n.tax+1):(2*n.tax-1))==min.node), which(((n.tax+1):(2*n.tax-1))==med.node)))] , num.cal )

        }
    }
# END OF ADDITIONAL SECTION



		cal.times <- matrix(NA, nrow=length(cal.nodes), ncol=1)
		rownames(cal.times) <- cal.nodes

		cal.tax <- list()

		for(a in 1:length(cal.nodes)){
			cal.times[a] <- br.times.m[rownames(br.times.m)==cal.nodes[a]]
			cal.tax[[a]] <- tips(sim.dat[[2]], cal.nodes[a])
		}
		names(cal.tax) <- cal.nodes
		#saving to calib.sim
		calib.sim <- list()
	#calib.sim[[1]] is the calibrations
		calib.sim[[1]] <- list()
		calib.sim[[1]][[1]] <- num.cal
		calib.sim[[1]][[2]] <- cal.times
		calib.sim[[1]][[3]] <- list()
		for(b in 1:length(cal.tax)){
			calib.sim[[1]][[3]][[b]] <- cal.tax[[b]]
		}
	#calib.sim[[2]] is the minimum and median nodes
		calib.sim[[2]] <- list()
		calib.sim[[2]][[1]] <- "min and median nodes"
		calib.sim[[2]][[2]] <- matrix(c(min.time, med.time), 2, 1)
		rownames(calib.sim[[2]][[2]]) <- c(min.node, med.node)
		calib.sim[[2]][[3]] <- list()
		for(c in 1:2){
			calib.sim[[2]][[3]][[c]] <- interest.tax[[c]]
		}

	if(calib.sim[[1]][[1]]==1 && calib.sim[[1]][[2]]==max(br.times.m)){
		calib.sim[[3]] <- "SingleCalibrationAtRoot"
		}else if(calib.sim[[1]][[1]]==1 && calib.sim[[1]][[2]] != max(br.times.m)){
			calib.sim[[3]] <- "SingleCalibrationNotAtRoot"
			}else if(calib.sim[[1]][[1]] > 1 && any(calib.sim[[1]][[2]]==max(br.times.m))){
				calib.sim[[3]] <- "MultipleCalibrationsOneAtRoot"
				}else if(calib.sim[[1]][[1]] >1 && !any(calib.sim[[1]][[2]]==max(br.times.m))){
					calib.sim[[3]] <- "MultipleCalibrationsNoneAtRoot"
					}

	# Print message of weather the mean or median are nested within any of the calibrations

	min.in.calib <- vector()
	med.in.calib <- vector()
	for(c in 1:num.cal){
		min.in.calib[c] <- all(min.tax %in% calib.sim[[1]][[3]][[c]])
		med.in.calib[c] <- all(med.tax %in% calib.sim[[1]][[3]][[c]])
	}

	if(any(min.in.calib) && any(med.in.calib)){
		calib.sim[[4]] <- "min&medNodesNested"
	}else if(any(min.in.calib) && !any(med.in.calib)){
			calib.sim[[4]] <- "minNodeNested"
			}else if(!any(min.in.calib) && any(med.in.calib)){
					calib.sim[[4]] <- "medNodeNested"
					}else if(!any(min.in.calib) && !any(med.in.calib)){
							calib.sim[[4]] <- "NoNodesNested"
							}


	names(calib.sim) <- c("calib_info", "target_node_info", "cals_at_root", "nestedness")
	
	names(calib.sim[[1]]) <- c("n_cals", "age_cals", "nodes_cals")
	
	
	return(calib.sim)


}
