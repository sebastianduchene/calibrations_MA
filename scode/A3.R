# Version edited at 11h34 of Wed March 20 2013
# 	Now it also reports the 95% HPD for the coefficient of rate variation.


run.beast.analysis <- function(paths){
	paths <- readLines(paths)
	###Run Beast

	if(!any(grep("xml", dir()))){
		print("No xml files here to run")
		}else {
			xml.file <- dir()[grep(".xml", dir())][1]

			if(any(grep(".log", dir()))){
				logs.rm <- grep(".log", dir(), value = T)
				for(a in 1:length(logs.rm)){
					system(paste("rm", logs.rm[a]))
					}
				}

			system(paste(paths[1], xml.file))

			# BEAST STARTS RUNNING HERE

			system(paste(paths[2], "simulated.data.log", "analyzed.txt"))

# Read analyses from analyzed.txt and process into a data matrix

###
	lineas <- readLines("analyzed.txt")
	if(length(lineas!=0)){
		if(any(grepl("distributionIndex", lineas))){
			ending <- grep("distributionIndex", lineas)-2
			}else{
				ending <- grep("coefficientOf", lineas)-1
			}
	res <- read.table("analyzed.txt", skip=1, head=T, as.is=T, nrow=ending, comment.char="*")
###
		}
	res.matrix.names <- c("posterior", "posteriorESS", "prior","priorESS", "likelihood", 	"likelihoodESS",

"rootHeight", "rootHeight95Lower", "rootHeight95Upper", "rootHeight50Lower", "rootHeight50Upper", "rootCalibrated", "RealRootTMRCA","rootHeightESS",
 "tmrca(cal1)","tmrca(cal1)95Lower","tmrca(cal)95Higher","tmrca(cal1)50Lower","tmrca(cal)50Higer","tmrca(cal1)Prior","tmrca(cal1)ESS",

"tmrca(min)","tmrca(min)95Lower","tmrca(min)95Higher", "tmrca(min)50Lower","tmrca(min)50Higher","RealTMRCA(min)","tmrca(min)ESS",

"tmrca(med)","tmrca(med)95Lower","tmrca(med)95Higher", "tmrca(med)50Lower","tmrca(med)50Higher","RealTMRCA(med)","tmrca(med)ESS",

"meanRate","meanRate95Low", "meanRate95High", "meanRate50Low", "meanRate50High", "realMeanRate","meanRateESS",

"CoeffOfVar", "CoeffOfVar95High", "CoeffOfVar95Low","CoeffOfVarESS",

"Covar", "CovarESS", "min.nested", "med.nested", "clock.mean", "clock.mean50Low","clock.mean50High","clock.mean95Low","clock.mean95High","clock.sd", "log.clock", "distributionIndex")


	res.matrix <- matrix(NA, nrow = 1, ncol = 58)
	colnames(res.matrix) <- res.matrix.names
	if(length(lineas)==0){
		res.matrix[1,1:50] <- NA
	}else{

	#Prior, posterior etc
	res.matrix[1,1] <- res$mean[1]
	res.matrix[1,2] <- res$ESS[1]
	res.matrix[1,3] <- res$mean[2]
	res.matrix[1,4] <- res$ESS[2]
	res.matrix[1,5] <- res$mean[3]
	res.matrix[1,6] <- res$ESS[3]
	#tmrcaRoot
	res.matrix[1,7] <- res[4,2]
	res.matrix[1,8] <- res[4,5]
	res.matrix[1,9] <- res[4,6]
	res.matrix[1,10] <- res[4,8]
	res.matrix[1,11] <- res[4,9]
	if(calib.sim[[3]] == "SingleCalibrationAtRoot"){
		res.matrix[1,12] <- 1
		}else{
			res.matrix[1,12] <- 0
			}
	res.matrix[1,13] <- max(branching.times(sim.dat[[2]]))
	res.matrix[1,14] <- res[4,7]


	if(length(calib.sim[[1]][[3]]) == 1){
	#tmrca cal
		res.matrix[1,15] <- res[5,2]
		res.matrix[1,16] <- res[5,5]
		res.matrix[1,17] <- res[5,6]
		res.matrix[1,18] <- res[5,8]
		res.matrix[1,19] <- res[5,9]
		res.matrix[1,20] <- calib.sim[[1]][[2]][1]
		res.matrix[1,21] <- res[5,7]
		}else if(length(calib.sim[[1]][[3]]) > 1){
			res.matrix[1,15] <- median(calib.sim[[1]][[2]])
			res.matrix[1,16] <- length(calib.sim[[1]][[3]])
			colnames(res.matrix)[15:21] <- c("medianCalDepth", "Num.Calibrations", "N", "N", "N", "N", "N")
			}
#tmrca min
	res.matrix[1,22] <- res[res[,1]=="tmrca(min)",2]
	res.matrix[1,23] <- res[res[,1]=="tmrca(min)",5]
	res.matrix[1,24] <- res[res[,1]=="tmrca(min)",6]
	res.matrix[1,25] <- res[res[,1]=="tmrca(min)",8]
	res.matrix[1,26] <- res[res[,1]=="tmrca(min)",9]
	res.matrix[1,27] <- calib.sim[[2]][[2]][1]
	res.matrix[1,28] <- res[res[,1]=="tmrca(min)",7]
	#tmrca med
	res.matrix[1,29] <- res[res[,1]=="tmrca(med)",2]
	res.matrix[1,30] <- res[res[,1]=="tmrca(med)",5]
	res.matrix[1,31] <- res[res[,1]=="tmrca(med)",6]
	res.matrix[1,32] <- res[res[,1]=="tmrca(med)",8]
	res.matrix[1,33] <- res[res[,1]=="tmrca(med)",9]
	res.matrix[1,34] <- calib.sim[[2]][[2]][2]
	res.matrix[1,35] <- res[res[,1]=="tmrca(med)",7]
	#rate
	res.matrix[1,36] <- res[res[,1]=="meanRate",2]
	res.matrix[1,37] <- res[res[,1]=="meanRate",5]
	res.matrix[1,38] <- res[res[,1]=="meanRate",6]
	res.matrix[1,39] <- res[res[,1]=="meanRate",8]
	res.matrix[1,40] <- res[res[,1]=="meanRate",9]
	res.matrix[1,41] <- exp(1)^sim.dat[[4]]$rate
	res.matrix[1,42] <- res[res[,1]=="meanRate",7]
	# Coeffs
	res.matrix[1,43] <- res[res[,1]=="coefficientOfVariation",2]
	res.matrix[1,44] <- res[res[,1]=="coefficientOfVariation",6]
	res.matrix[1,45] <- res[res[,1]=="coefficientOfVariation",5]

	res.matrix[1,46] <- res[res[,1]=="coefficientOfVariation",7]
	res.matrix[1,47] <- res[res[,1]=="covariance",2]
	res.matrix[1,48] <- res[res[,1]=="covariance",7]
	#nodes of interest nested
	if(calib.sim[[4]] == "min&medNodesNested"){
		res.matrix[1,49] <- 1
		res.matrix[1,50] <- 1
		}else if(calib.sim[[4]] == "minNodeNested"){
			res.matrix[1,49] <- 1
			res.matrix[1, 50] <- 0
			}else if(calib.sim[[4]] == "medNodeNested"){
				res.matrix[1,49] <- 0
				res.matrix[1,50] <- 1
				}else{
					res.matrix[1,49] <- 0
					res.matrix[1,50] <- 0
					}
		res.matrix[1,51] <- res[which(res[,1]=="ucld.mean" | res[,1]=="uced.mean")[1] , 2]
		res.matrix[1,52] <- res[which(res[,1]=="ucld.mean" | res[,1]=="uced.mean")[1] , 8]
		res.matrix[1,53] <- res[which(res[,1]=="ucld.mean" | res[,1]=="uced.mean")[1] , 9]
		res.matrix[1,54] <- res[which(res[,1]=="ucld.mean" | res[,1]=="uced.mean")[1] , 5]
		res.matrix[1,55] <- res[which(res[,1]=="ucld.mean" | res[,1]=="uced.mean")[1] , 6]
		res.matrix[1,56] <- res[which(res[,1]=="ucld.stdev" | res[,1]=="uced.mean")[1] , 2]

	if(any(res[,1]=="uced.mean")){
		colnames(res.matrix)[56] <- "uced.mean"
		}

		if(sim.dat[[4]]$sdrate==0){
			res.matrix[1,57] <- 0
			}else{
				res.matrix[1,57] <- 1
				}
	if(any(grepl("distributionIndex", lineas))){
            print("This is a mixture model analysis")
            #system("sleep 3")
		res.matrix[1, 58] <- res[res[,1]=="branchRates.distributionIndex", 2]
		}else{
                    print("This is not a mixture model analysis")
             #       system("sleep 3")
                }


}
	if(any(dir()=="results.matrix.csv")){
		write.table(res.matrix, "results.matrix.csv", row.names=F, col.names=F, append=T,sep=",")
		}else{
			write.table(res.matrix, "results.matrix.csv", row.names=F, append=F, sep=",")
			}


		}

		system("rm simulated.data.log analyzed.txt")
		system(paste("rm", xml.file))

}





