



###########################
###########################
###########################
###########################
###########################
#source("~/Desktop/CalibProjectScripts/scode/A4.EXECUTE.R")
#run.simulations(1000, 20, -2.3, 0.1, 50, F, F, 1000, 100, "paths.txt", 10)

#run.simulations <- function(seq.length, n.tax, mean.rate, sd.rate, root.height, mult.cal, exclude.root, mcmc.length, mcmc.freq, paths, n.simulations)

#run.simulations <- function(seq.length, n.tax, mean.rate, sd.rate, root.height, mult.cal, #exclude.root, mcmc.length, mcmc.freq, paths, n.simulations){
###########################
###########################
###########################
###########################
###########################





# Settings

#	seq.length = 2000		#length of sequences to be simulated
#	n.tax = 50 				#number of taxa to simulate			
#	mean.rate = -6.908		#rate of substitutions for sequence simulation
#	sd.rate = 0				#standard deviation of the rate in the simulations. Set to 0 for a strict clock
#	root.height = 50		#Height of the root. This is the TMRCA of all the taxa in the tree
#	mult.cal = T			#Select T to include multiple calibrations, or F for single calibration, and "fixed" to fix the number of calibrations
#	n.cal = 5				#Put the number of calibrations for the multiple calibration analysis.
#	exclude.root = F		#Select T if the root should be excluded as a possible calibration, select F otherwise
#	rate.dist = "lognormal"	#distribution for the simulated rate, can be lognormal or exponential. For strict clock select "lognorma" and set sd.rate=0
#	clock = "lognormal"		#clock model for analysis. The options are lognormal or exponential
#	mcmc.length = 2000000	#length of the run
#	mcmc.freq = 1000		#log frequency. This value applies to the log on the screen, and the logfile
#	paths="paths.txt"		#name of the file where the paths for beast and loganalyzer binaries are found
#	n.simulations = 1000	#number of simulations to run. Each time a new dataset is generated and analyzed

#############################
#############################
#############################
#############################
#############################
#############################
	
#	setwd("..")
	setwd("../scode")
	source("A0.R")
	source("A1.R")
	source("A2.R")
	source("A3.R")
	setwd("..")
	setwd(settings$settings[i])
	
	a <- 0
	while(a <= n.simulations){
			
		#Removing any old files
		if(any(grep(".xml", dir()))){
			l <- length((grep(".xml", dir())))
			for(b in 1:l){
				system(paste("rm", dir()[grep(".xml", dir())][1]))
			}
		}
		if(any(grep(".log", dir()))){
			l <- length((grep(".log", dir())))
			for(b in 1:l){
				system(paste("rm", dir()[grep(".log", dir())][1]))
			}
		}
		if(any(grep("analyzed.txt", dir()))){
			l <- length((grep("analyzed.txt", dir())))
			for(b in 1:l){
				system(paste("rm", dir()[grep("analyzed.txt", dir())][1]))
			}
		}
	sim.dat <- sim.dataset(seq.length, n.tax, mean.rate, sd.rate, rate.dist ,root.height)
	calib.sim <- select.calibration(sim.dat, mult.cal, n.cal, exclude.root)
	generate.xml.file(sim.dat, calib.sim, clock ,mcmc.length, mcmc.freq, mix.model)
	run.beast.analysis(paths)
	
	remove(sim.dat, calib.sim)
	
	print(paste("ANALYSIS", a+1 , "OF", n.simulations, ";" , 100*((a+1)/n.simulations), "%"))

	a <- a+1	

}
	
