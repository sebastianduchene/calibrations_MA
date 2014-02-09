# Modified in 11h37 Wed March 20 2013 needed to nest the settings so that it doesn't halt aftere running one of them

# First version for multiple calibrations with fixed number. 00h39 Tuesday Mar 19 2013
                                        # This is the code to rdun several of the calibration settings.

# First we set the working directory to one level higher and read the file with the settings

settings.all <- read.table("settings.csv", sep = ",", head = T, as.is = T)


# Now define how which of the settings are to be run in this analysis

#settings <- settings[1:4,]

settings.subset <- settings.all[c(9:15),]


# To run those available

#settings.done <- grep("settings*.", dir(), value = T)

#settings.remaining <- settings$settings[!(settings$settings %in% settings.done)]

#settings <- settings[settings$settings %in% settings.remaining, ]
#for(i in 1:length(settings.all)){

#	settings.done <- grep("settings*.", dir(), value = T)

#	settings.remaining <- settings$settings[!(settings$settings %in% settings.done)]

#	settings <- settings[settings$settings %in% settings.remaining, ]

#for(i in 1:nrow(settings.subset)){

while(length(grep("settings*.", dir(), value = T)) < nrow(settings.subset)){

	settings.done <- grep("settings*.", dir(), value = T)
	settings.remaining <- settings.subset$settings[!(settings.subset$settings %in% settings.done)]

#	if(length(settings.remaninig) != length(settings.done)){
		settings <- settings.subset[settings.subset$settings %in% settings.remaining, ]

	i=1

    	system(paste("mkdir", settings$settings[i]))
    	system(paste("cp", "./scode/A4.set.runs.R", settings$settings[i]))
    	system(paste("cp", paste("./", settings$paths[i], sep =""), paste(settings$settings[i], "/", settings$paths, sep = "" )))
    	setwd(settings$settings[i])

		seq.length = settings$seq.length[i]		#length of sequences to be simulated
		n.tax = settings$n.tax[i] 			#number of taxa to simulate
        mean.rate = settings$mean.rate[i]		#rate of substitutions for sequence simulation
		sd.rate = settings$sd.rate[i]			#standard deviation of the rate in the simulations. Set to 0 for a strict clock
        root.height = settings$root.height[i]		#Height of the root. This is the TMRCA of all the taxa in the tree
		mult.cal = settings$mult.cal[i]			#Select T to include multiple calibrations, F for single calibration, or "fixed" to fix the number of calibrations
		n.cal = settings$n.cal[i]			#Put the number of calibrations for the multiple calibration analysis.
        exclude.root = settings$exclude.root[i]         #Select T if the root should be excluded as a possible calibration, select F otherwise
		rate.dist = settings$rate.dist[i]	        #distribution for the simulated rate, can be lognormal or exponential. For strict clock select "lognormal" and set sd.rate=0
		clock = settings$clock[i]		        #clock model for analysis. The options are lognormal or exponential
		mcmc.length = settings$mcmc.length[i]	        #length of the run
		mcmc.freq = settings$mcmc.freq[i]		#log frequency. This value applies to the log on the screen, and the logfile
		paths = settings$paths[i]	 	        #name of the file where the paths for beast and loganalyzer binaries are found
		n.simulations = settings$n.simulations[i]       #Number of simulations to run for this analysis
		mix.model = settings$mix.model[i]	#Select T to run model averaging

		source("../scode/A4.set.runs.R")
	
		system("rm A4.set.runs.R")

    	setwd("..")
#		}else{
#			print("No more simulation settings to run")
#			}
	

}

print(" Analysis for all settings complete")


