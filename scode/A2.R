# Calibration and tmrca estimation project script A.2
# Modified Feb 8 2014- The xml file can also include mixture models for lognormal and exponential clock models
# Modified September 24 2012- Now the xml file can be written to use beast's implementations of lognormal or exponential clocks


generate.xml.file <- function(sim.dat, calib.sim, clock, mcmc.length, mcmc.freq, mix.model = F){
	#General variables for later use
	tree <- sim.dat[[2]]
	secs <- sim.dat[[3]]
	n.cals <- calib.sim[[1]][[1]]
	tax.cals <- calib.sim[[1]][[3]]
	cal.names <- paste("cal", 1:n.cals, sep="")
	min.med.names <- c("min", "med")
	tax.min.med <- calib.sim[[2]][[3]]
	cal.times <- calib.sim[[1]][[2]]
	#Begin with the first lines required for beast
	
#	block.1 <- scan(,what="")
	block.1 <- character()
	block.1[length(block.1)+1]<-"<?xml version=\"1.0\" standalone=\"yes\"?>"
	block.1[length(block.1)+1]<-"<!-- Generated for BEAST 1.7.2- 2012                                        -->"
	block.1[length(block.1)+1]<-"<!--                -->"
	block.1[length(block.1)+1]<-"<!--               -->"
	block.1[length(block.1)+1]<-"<!--       I        -->"
	block.1[length(block.1)+1]<-"<!--       -->"
	block.1[length(block.1)+1]<-"<!--                                             -->"
	block.1[length(block.1)+1]<-"<beast>"
	block.1[length(block.1)+1]<-"<taxa id=\"all.taxa\">"

# Now we write the block for all the taxa
	tax.begin <- "<taxon id=\""
	tax.end <- "\"/>"
	block.taxa <- character()
	for(a in 1:length(row.names(secs))){
		block.taxa[a] <- paste(tax.begin, row.names(secs)[a], tax.end, sep="")
	}
	block.taxa[length(block.taxa)+1] <- "</taxa>"

#Now the blocks for each calibration	
	block.calib.tax <- character()
	for(b in 1:n.cals){
		block.calib.tax[length(block.calib.tax)+1]  <- paste("<taxa id=\"",cal.names[b], "\">" , sep="")
		for(c in 1:length(tax.cals[[b]])){
			block.calib.tax[length(block.calib.tax)+1] <- paste("<taxon idref=\"",tax.cals[[b]][c], tax.end, sep="")
		}
		block.calib.tax[length(block.calib.tax)+1] <- "</taxa>"		
	}
	
#Now the blocks for the min and median nodes

	block.min.med <- character()	
	for(d in 1:2){
		block.min.med[length(block.min.med)+1] <- paste("<taxa id=\"", min.med.names[d], "\">", sep="")
		for(e in 1:length(tax.min.med[[d]])){
			block.min.med[length(block.min.med)+1] <- paste("<taxon idref=\"", tax.min.med[[d]][e], tax.end, sep="")
		}
		block.min.med[length(block.min.med)+1] <- "</taxa>"			
	}

# Now we need the block for the sequences

	block.secs <- character()
	block.secs[length(block.secs)+1] <- "	<alignment id=\"alignment\" dataType=\"nucleotide\">"

	for(f in 1:nrow(secs)){
		block.secs[length(block.secs)+1] <- "<sequence>"
		block.secs[length(block.secs)+1] <- paste("<taxon idref=\"", rownames(secs)[f], "\"/>", sep="")
		block.secs[length(block.secs)+1] <- paste(as.character(secs[f,]), collapse="", sep="")
		block.secs[length(block.secs)+1] <- "</sequence>"
	}
	block.secs[length(block.secs)+1] <- "</alignment>"

	block.pat.tree.model.size <- character()
	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "	<patterns id=\"patterns\" from=\"1\" strip=\"false\">"
	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "		<alignment idref=\"alignment\"/>"
	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "	</patterns>"
	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "	<yuleModel id=\"yule\" units=\"substitutions\">"
	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "		<birthRate>"
	if(clock=="lognormal"){
	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- paste("<parameter id=\"yule.birthRate\" value=\"", exp(1)^sim.dat[[4]]$rate ,  "\" lower=\"0.0\"/>", sep="")
	}else if(clock=="exponential"){
			block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- paste("<parameter id=\"yule.birthRate\" value=\"", 0.3333 ,  "\" lower=\"0.0\"/>", sep="")
			}
	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "		</birthRate>"
	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "	</yuleModel>"
#	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "	<constantSize id=\"initialDemo\" units=\"substitutions\">"
#	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "		<populationSize>"
#	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "			<parameter id=\"initialDemo.popSize\" value=\"100.0\"/>"
#	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "		</populationSize>"
#	block.pat.tree.model.size[length(block.pat.tree.model.size)+1] <- "	</constantSize>"

	block.tree.fix <- character()
	block.tree.fix[length(block.tree.fix)+1] <- "<newick id=\"startingTree\" usingDates=\"false\">"
	block.tree.fix[length(block.tree.fix)+1] <- write.tree(tree)
	block.tree.fix[length(block.tree.fix)+1] <- "</newick>"

	block.tree.param <- character()
	block.tree.param[length(block.tree.param)+1] <-"<treeModel id=\"treeModel\">"
	block.tree.param[length(block.tree.param)+1] <-"<tree idref=\"startingTree\"/>"
	block.tree.param[length(block.tree.param)+1] <-"<rootHeight>"
	block.tree.param[length(block.tree.param)+1] <-"<parameter id=\"treeModel.rootHeight\"/>"
	block.tree.param[length(block.tree.param)+1] <-"</rootHeight>"
	block.tree.param[length(block.tree.param)+1] <-"<nodeHeights internalNodes=\"true\">"
	block.tree.param[length(block.tree.param)+1] <-"<parameter id=\"treeModel.internalNodeHeights\"/>"
	block.tree.param[length(block.tree.param)+1] <-"</nodeHeights>"
	block.tree.param[length(block.tree.param)+1] <-"<nodeHeights internalNodes=\"true\" rootNode=\"true\">"
	block.tree.param[length(block.tree.param)+1] <-"<parameter id=\"treeModel.allInternalNodeHeights\"/>"
	block.tree.param[length(block.tree.param)+1] <-"</nodeHeights>"
	block.tree.param[length(block.tree.param)+1] <-"</treeModel>"


#Now taxon sets statistics for calibrations
	block.cal.stat <- character()
	for(g in 1:n.cals){
		block.cal.stat[length(block.cal.stat)+1] <- paste("<tmrcaStatistic id=\"tmrca(",cal.names[g],")\" includeStem=\"false\">", sep="")
		block.cal.stat[length(block.cal.stat)+1] <- "<mrca>"
		block.cal.stat[length(block.cal.stat)+1] <- paste("<taxa idref=\"", cal.names[g],"\"/>", sep="")
		block.cal.stat[length(block.cal.stat)+1] <- "</mrca>"
		block.cal.stat[length(block.cal.stat)+1] <-"<treeModel idref=\"treeModel\"/>"
		block.cal.stat[length(block.cal.stat)+1] <-"</tmrcaStatistic>"		
	}
	
	block.min.med.stat <- character()
	for(h in 1:2){
		block.min.med.stat[length(block.min.med.stat)+1] <- paste("<tmrcaStatistic id=\"tmrca(",min.med.names[h],")\" includeStem=\"false\">", sep="")
		block.min.med.stat[length(block.min.med.stat)+1] <- "<mrca>"
		block.min.med.stat[length(block.min.med.stat)+1] <- paste("<taxa idref=\"", min.med.names[h],"\"/>", sep="")
		block.min.med.stat[length(block.min.med.stat)+1] <- "</mrca>"
		block.min.med.stat[length(block.min.med.stat)+1] <-"<treeModel idref=\"treeModel\"/>"
		block.min.med.stat[length(block.min.med.stat)+1] <-"</tmrcaStatistic>"
	}
	
	block.spec <- character()
	block.spec[length(block.spec)+1] <- "<speciationLikelihood id=\"speciation\">"
	block.spec[length(block.spec)+1] <- "<model>"
	block.spec[length(block.spec)+1] <- "<yuleModel idref=\"yule\"/>"
	block.spec[length(block.spec)+1] <-	"</model>"
	block.spec[length(block.spec)+1] <- "<speciesTree>"
	block.spec[length(block.spec)+1] <-	"<treeModel idref=\"treeModel\"/>"
	block.spec[length(block.spec)+1] <- "</speciesTree>"
	block.spec[length(block.spec)+1] <- "</speciationLikelihood>"


#####	
	block.clock <- character()
	if(mix.model == T){
		block.clock[length(block.clock)+1] <- "<mixtureModelBranchRates id=\"branchRates\">"
		block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
		block.clock[length(block.clock)+1] <-	"<!-- INSERTED FIRST CLOCK -->"
		block.clock[length(block.clock)+1] <- "<distribution>"
		block.clock[length(block.clock)+1] <- "<logNormalDistributionModel meanInRealSpace=\"true\" stdevInRealSpace=\"true\">"
		block.clock[length(block.clock)+1] <- "<mean>"
		block.clock[length(block.clock)+1] <- "<parameter id=\"ucld.mean\" value=\"1.0\"/>"
		block.clock[length(block.clock)+1] <- "</mean>"
		block.clock[length(block.clock)+1] <- "<stdev>"
		block.clock[length(block.clock)+1] <- "<parameter id=\"ucld.stdev\" value=\"0.1\" lower=\"0.0\" upper=\"10.0\"/>"
		block.clock[length(block.clock)+1] <- "</stdev>"
		block.clock[length(block.clock)+1] <- "</logNormalDistributionModel>"
		block.clock[length(block.clock)+1] <- "</distribution>"
		block.clock[length(block.clock)+1] <- "<!--INSERT SECOND CLOCK HERE  -->" 
		block.clock[length(block.clock)+1] <- "<distribution>"
		block.clock[length(block.clock)+1] <- "<exponentialDistributionModel>"
		block.clock[length(block.clock)+1] <- "<mean>"
		block.clock[length(block.clock)+1] <- "<parameter id=\"uced.mean\" value=\"1.0\"/>"
		block.clock[length(block.clock)+1] <- "</mean>"
		block.clock[length(block.clock)+1] <- "</exponentialDistributionModel>"
		block.clock[length(block.clock)+1] <- "</distribution>"
		block.clock[length(block.clock)+1] <- "<distributionIndex>"
		block.clock[length(block.clock)+1] <- "<parameter id=\"branchRates.distributionIndex\"/>"
		block.clock[length(block.clock)+1] <- "</distributionIndex>"
		block.clock[length(block.clock)+1] <- "<rateCategoryQuantiles>"
		block.clock[length(block.clock)+1] <- "<parameter id=\"branchRates.categoryQuantiles\" dimension=\"22\" value=\"0.5\" lower=\"0.01\" upper=\"0.99\"/>"
		block.clock[length(block.clock)+1] <- "</rateCategoryQuantiles>"
		block.clock[length(block.clock)+1] <- "</mixtureModelBranchRates>"

		block.clock[length(block.clock)+1] <- "<rateStatistic id=\"meanRate\" name=\"meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">"
		block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
		block.clock[length(block.clock)+1] <- "<mixtureModelBranchRates idref=\"branchRates\"/>"
		block.clock[length(block.clock)+1] <- "</rateStatistic>"

		block.clock[length(block.clock)+1] <- "<rateStatistic id=\"coefficientOfVariation\" name=\"coefficientOfVariation\" mode=\"coefficientOfVariation\" internal=\"true\" external=\"true\">"
		block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
		block.clock[length(block.clock)+1] <- "<!-- REPLACE DISCRETIZED BRANCH RATES BU MIXTURE MODEL BRANCH RATES -->"
		block.clock[length(block.clock)+1] <- "<mixtureModelBranchRates idref=\"branchRates\"/>"
		block.clock[length(block.clock)+1] <- "</rateStatistic>"

		block.clock[length(block.clock)+1] <- "<rateCovarianceStatistic id=\"covariance\" name=\"covariance\">"
		block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
		block.clock[length(block.clock)+1] <- "<mixtureModelBranchRates idref=\"branchRates\"/>"
		block.clock[length(block.clock)+1] <- "</rateCovarianceStatistic>"
	}


	if(clock=="lognormal" && mix.model==F){
		block.clock[length(block.clock)+1] <- "<discretizedBranchRates id=\"branchRates\">"
		block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
		block.clock[length(block.clock)+1] <- "<distribution>"
		block.clock[length(block.clock)+1] <- "<logNormalDistributionModel meanInRealSpace=\"true\">"
		block.clock[length(block.clock)+1] <- "<mean>"
		block.clock[length(block.clock)+1] <- paste("<parameter id=\"ucld.mean\" value=\"", exp(1)^sim.dat[[4]]$rate ,  "\" lower=\"0.0\"/>", sep="")
		block.clock[length(block.clock)+1] <- "</mean>"
		block.clock[length(block.clock)+1] <- "<stdev>"
		block.clock[length(block.clock)+1] <- "<parameter id=\"ucld.stdev\" value=\"0.3333333333333333\" lower=\"0.0\"/>"
		block.clock[length(block.clock)+1] <- "</stdev>"
		block.clock[length(block.clock)+1] <- "</logNormalDistributionModel>"
		block.clock[length(block.clock)+1] <- "</distribution>"
		block.clock[length(block.clock)+1] <- "<rateCategories>"
		block.clock[length(block.clock)+1] <- "<parameter id=\"branchRates.categories\"/>"
		block.clock[length(block.clock)+1] <- "</rateCategories>"
		block.clock[length(block.clock)+1] <- "</discretizedBranchRates>"
		block.clock[length(block.clock)+1] <- "<rateStatistic id=\"meanRate\" name=\"meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">"
		block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
		block.clock[length(block.clock)+1] <- "<discretizedBranchRates idref=\"branchRates\"/>"
		block.clock[length(block.clock)+1] <- "</rateStatistic>"
		block.clock[length(block.clock)+1] <- "<rateStatistic id=\"coefficientOfVariation\" name=\"coefficientOfVariation\" mode=\"coefficientOfVariation\" internal=\"true\" external=\"true\">"
		block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
		block.clock[length(block.clock)+1] <- "<discretizedBranchRates idref=\"branchRates\"/>"
		block.clock[length(block.clock)+1] <- "</rateStatistic>"
		block.clock[length(block.clock)+1] <- "<rateCovarianceStatistic id=\"covariance\" name=\"covariance\">"
		block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
		block.clock[length(block.clock)+1] <- "<discretizedBranchRates idref=\"branchRates\"/>"
		block.clock[length(block.clock)+1] <- "</rateCovarianceStatistic>"
		}else if(clock=="exponential" && mix.model==F){
					block.clock[length(block.clock)+1] <- "<discretizedBranchRates id=\"branchRates\">"
					block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
					block.clock[length(block.clock)+1] <- "<distribution>"
					block.clock[length(block.clock)+1] <- "<exponentialDistributionModel>"
					block.clock[length(block.clock)+1] <- "<mean>"
					block.clock[length(block.clock)+1] <- paste("<parameter id=\"uced.mean\" value=\"",  1.0,  "\" lower=\"0.0\"/>", sep="")
					block.clock[length(block.clock)+1] <- "</mean>"
					block.clock[length(block.clock)+1] <- "</exponentialDistributionModel>"
					block.clock[length(block.clock)+1] <- "</distribution>"
					block.clock[length(block.clock)+1] <- "<rateCategories>"
					block.clock[length(block.clock)+1] <- "<parameter id=\"branchRates.categories\"/>"
					block.clock[length(block.clock)+1] <- "</rateCategories>"
					block.clock[length(block.clock)+1] <- "</discretizedBranchRates>"
					block.clock[length(block.clock)+1] <- "<rateStatistic id=\"meanRate\" name=\"meanRate\" mode=\"mean\" internal=\"true\" external=\"true\">"
					block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
					block.clock[length(block.clock)+1] <- "<discretizedBranchRates idref=\"branchRates\"/>"
					block.clock[length(block.clock)+1] <- "</rateStatistic>"
					block.clock[length(block.clock)+1] <- "<rateStatistic id=\"coefficientOfVariation\" name=\"coefficientOfVariation\" mode=\"coefficientOfVariation\" internal=\"true\" external=\"true\">"
					block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
					block.clock[length(block.clock)+1] <- "<discretizedBranchRates idref=\"branchRates\"/>"
					block.clock[length(block.clock)+1] <- "</rateStatistic>"
					block.clock[length(block.clock)+1] <- "<rateCovarianceStatistic id=\"covariance\" name=\"covariance\">"
					block.clock[length(block.clock)+1] <- "<treeModel idref=\"treeModel\"/>"
					block.clock[length(block.clock)+1] <- "<discretizedBranchRates idref=\"branchRates\"/>"
					block.clock[length(block.clock)+1] <- "</rateCovarianceStatistic>"
				}
		
### 
	block.model.ops <- character()
	
	block.model.ops[length(block.model.ops)+1] <-" <HKYModel id=\"hky\">"
	block.model.ops[length(block.model.ops)+1] <-"<frequencies>"
	block.model.ops[length(block.model.ops)+1] <-"<frequencyModel dataType=\"nucleotide\">"
	block.model.ops[length(block.model.ops)+1] <-"<frequencies>"
	block.model.ops[length(block.model.ops)+1] <-"<parameter id=\"frequencies\" value=\"0.25 0.25 0.25 0.25\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</frequencies>"
	block.model.ops[length(block.model.ops)+1] <-"</frequencyModel>"
	block.model.ops[length(block.model.ops)+1] <-"</frequencies>"
	block.model.ops[length(block.model.ops)+1] <-"<kappa>"
	block.model.ops[length(block.model.ops)+1] <-"<parameter id=\"kappa\" value=\"2.0\" lower=\"0.0\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</kappa>"
	block.model.ops[length(block.model.ops)+1] <-"</HKYModel>"
	block.model.ops[length(block.model.ops)+1] <-"<siteModel id=\"siteModel\">"
	block.model.ops[length(block.model.ops)+1] <-"<substitutionModel>"
	block.model.ops[length(block.model.ops)+1] <-"<HKYModel idref=\"hky\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</substitutionModel>"
	block.model.ops[length(block.model.ops)+1] <-"</siteModel>"
	block.model.ops[length(block.model.ops)+1] <-"<treeLikelihood id=\"treeLikelihood\" useAmbiguities=\"false\">"
	block.model.ops[length(block.model.ops)+1] <-"<patterns idref=\"patterns\"/>"
	block.model.ops[length(block.model.ops)+1] <-"<treeModel idref=\"treeModel\"/>"
	block.model.ops[length(block.model.ops)+1] <-"<siteModel idref=\"siteModel\"/>"
	block.model.ops[length(block.model.ops)+1] <-"<discretizedBranchRates idref=\"branchRates\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</treeLikelihood>"
	block.model.ops[length(block.model.ops)+1] <-"<operators id=\"operators\">"
	block.model.ops[length(block.model.ops)+1] <-"<scaleOperator scaleFactor=\"0.75\" weight=\"0.1\">"
	block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"kappa\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</scaleOperator>"
	block.model.ops[length(block.model.ops)+1] <-"<deltaExchange delta=\"0.01\" weight=\"0.1\">"
	block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"frequencies\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</deltaExchange>"
	if(clock=="lognormal" && mix.model==F){
		block.model.ops[length(block.model.ops)+1] <-"<scaleOperator scaleFactor=\"0.75\" weight=\"3\">"
		block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"ucld.mean\"/>"
		block.model.ops[length(block.model.ops)+1] <-"</scaleOperator>"
		block.model.ops[length(block.model.ops)+1] <-"<scaleOperator scaleFactor=\"0.75\" weight=\"3\">"
		block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"ucld.stdev\"/>"
		block.model.ops[length(block.model.ops)+1] <-"</scaleOperator>"
		}else if(clock=="exponential" && mix.model==F){
				block.model.ops[length(block.model.ops)+1] <- "<scaleOperator scaleFactor=\"0.75\" weight=\"3\">"
				block.model.ops[length(block.model.ops)+1] <- "<parameter idref=\"uced.mean\"/>"
				block.model.ops[length(block.model.ops)+1] <- "</scaleOperator>"
				
				}else if(mix.model==T){
					block.model.ops[length(block.model.ops)+1] <- "<scaleOperator scaleFactor=\"0.75\" weight=\"3\">"
					block.model.ops[length(block.model.ops)+1] <- "<parameter idref=\"ucld.stdev\"/>"
					block.model.ops[length(block.model.ops)+1] <- "</scaleOperator>"
				}
	block.model.ops[length(block.model.ops)+1] <-"<scaleOperator scaleFactor=\"0.75\" weight=\"3\">"
	block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"treeModel.rootHeight\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</scaleOperator>"
	block.model.ops[length(block.model.ops)+1] <-"<uniformOperator weight=\"30\">"
	block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"treeModel.internalNodeHeights\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</uniformOperator>"
	block.model.ops[length(block.model.ops)+1] <-"<scaleOperator scaleFactor=\"0.75\" weight=\"3\">"
	block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"yule.birthRate\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</scaleOperator>"
	block.model.ops[length(block.model.ops)+1] <-"<upDownOperator scaleFactor=\"0.75\" weight=\"3\">"
	
	if(clock=="lognormal" || mix.model==T){
		block.model.ops[length(block.model.ops)+1] <-"<up>"
		block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"ucld.mean\"/>"
		block.model.ops[length(block.model.ops)+1] <-"</up>"
		}else if(clock=="exponential" && mix.model == F){
					block.model.ops[length(block.model.ops)+1] <- "<up>"
					block.model.ops[length(block.model.ops)+1] <- "<parameter idref=\"uced.mean\"/>"
					block.model.ops[length(block.model.ops)+1] <- "</up>"
				}
				
			
	block.model.ops[length(block.model.ops)+1] <-"<down>"
	block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"treeModel.allInternalNodeHeights\"/>"
	block.model.ops[length(block.model.ops)+1] <-"</down>"
	block.model.ops[length(block.model.ops)+1] <-"</upDownOperator>"
	
	if(mix.model==F){
		block.model.ops[length(block.model.ops)+1] <-"<swapOperator size=\"1\" weight=\"10\" autoOptimize=\"false\">"
		block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"branchRates.categories\"/>"
		block.model.ops[length(block.model.ops)+1] <-"</swapOperator>"
		block.model.ops[length(block.model.ops)+1] <-"<uniformIntegerOperator weight=\"10\">"
		block.model.ops[length(block.model.ops)+1] <-"<parameter idref=\"branchRates.categories\"/>"
		block.model.ops[length(block.model.ops)+1] <-"</uniformIntegerOperator>"
		block.model.ops[length(block.model.ops)+1] <-"</operators>"
		}else if(mix.model==T){
				block.model.ops[length(block.model.ops)+1] <- "<uniformIntegerOperator weight=\"5.0\">"
				block.model.ops[length(block.model.ops)+1] <- "<parameter idref=\"branchRates.distributionIndex\"/>"
				block.model.ops[length(block.model.ops)+1] <- "</uniformIntegerOperator>"
				block.model.ops[length(block.model.ops)+1] <- "<uniformOperator weight=\"30.0\">"
				block.model.ops[length(block.model.ops)+1] <- "<parameter idref=\"branchRates.categoryQuantiles\"/>"
				block.model.ops[length(block.model.ops)+1] <- "</uniformOperator>"
				block.model.ops[length(block.model.ops)+1] <- "</operators>"
			}
			
	block.mcmc <- character()
	block.mcmc[length(block.mcmc)+1] <- paste("<mcmc id=\"mcmc\" chainLength=\"", as.integer(mcmc.length), "\" autoOptimize=\"true\">", sep="")
	block.mcmc[length(block.mcmc)+1] <- "<posterior id=\"posterior\">"
	block.mcmc[length(block.mcmc)+1] <- "<prior id=\"prior\">"
	
	for(i in 1:n.cals){
		block.mcmc[length(block.mcmc)+1] <- paste("<normalPrior mean=\"", round(cal.times[i], 2), "\" stdev=\"", round((0.1*cal.times[i]), 2), "\">", sep="")
		block.mcmc[length(block.mcmc)+1] <- paste("<statistic idref=\"tmrca(", cal.names[i], ")\"/>", sep="")
		block.mcmc[length(block.mcmc)+1] <- "</normalPrior>"
	}
	
	
	block.mcmc[length(block.mcmc)+1] <-"<logNormalPrior mean=\"1.0\" stdev=\"1.25\" offset=\"0.0\" meanInRealSpace=\"false\">"
	block.mcmc[length(block.mcmc)+1] <-"<parameter idref=\"kappa\"/>"
	block.mcmc[length(block.mcmc)+1] <-"</logNormalPrior>"
	block.mcmc[length(block.mcmc)+1] <-"<uniformPrior lower=\"0.0\" upper=\"1.0\">"
	block.mcmc[length(block.mcmc)+1] <-	"<parameter idref=\"frequencies\"/>"
	block.mcmc[length(block.mcmc)+1] <-	"</uniformPrior>"
	
	#clock-specific statistics in mcmc
	if(clock=="lognormal" || mix.model==T){
		block.mcmc[length(block.mcmc)+1] <-	"<exponentialPrior mean=\"0.3333333333333333\" offset=\"0.0\">"
		block.mcmc[length(block.mcmc)+1] <-	"<parameter idref=\"ucld.stdev\"/>"
		block.mcmc[length(block.mcmc)+1] <-	"</exponentialPrior>"	
		block.mcmc[length(block.mcmc)+1] <-	"<uniformPrior lower=\"0.0\" upper=\"1.0E100\">"
		block.mcmc[length(block.mcmc)+1] <-	"<parameter idref=\"ucld.mean\"/>"
		block.mcmc[length(block.mcmc)+1] <-	"</uniformPrior>"
		}else if(clock=="exponential" && mix.model==F){
				block.mcmc[length(block.mcmc)+1] <- "<exponentialPrior mean=\"1.0\" offset=\"0.0\">"
				block.mcmc[length(block.mcmc)+1] <-	"<parameter idref=\"uced.mean\"/>"
				block.mcmc[length(block.mcmc)+1] <- "</exponentialPrior>"
				}

	block.mcmc[length(block.mcmc)+1] <-	"<uniformPrior lower=\"0.0\" upper=\"1.0E100\">"
	block.mcmc[length(block.mcmc)+1] <-	"<parameter idref=\"yule.birthRate\"/>"
	block.mcmc[length(block.mcmc)+1] <-	"</uniformPrior>"
	block.mcmc[length(block.mcmc)+1] <-	"<speciationLikelihood idref=\"speciation\"/>"
	block.mcmc[length(block.mcmc)+1] <-	"</prior>"
	block.mcmc[length(block.mcmc)+1] <-	"<likelihood id=\"likelihood\">"
	block.mcmc[length(block.mcmc)+1] <-	"<treeLikelihood idref=\"treeLikelihood\"/>"
	block.mcmc[length(block.mcmc)+1] <-	"</likelihood>"
	block.mcmc[length(block.mcmc)+1] <-	"</posterior>"
	block.mcmc[length(block.mcmc)+1] <-	"<operators idref=\"operators\"/>"
		
	block.screen.log <- character()
	block.screen.log[length(block.screen.log)+1] <- paste("<log id=\"screenLog\" logEvery=\"", as.integer(mcmc.freq) ,"\">", sep="")
	block.screen.log[length(block.screen.log)+1] <- "<column label=\"Posterior\" dp=\"4\" width=\"12\">"
	block.screen.log[length(block.screen.log)+1] <- "<posterior idref=\"posterior\"/>"
	block.screen.log[length(block.screen.log)+1] <- "</column>"
	block.screen.log[length(block.screen.log)+1] <- "<column label=\"Prior\" dp=\"4\" width=\"12\">"
	block.screen.log[length(block.screen.log)+1] <- "<prior idref=\"prior\"/>"
	block.screen.log[length(block.screen.log)+1] <- "</column>"
	block.screen.log[length(block.screen.log)+1] <- "<column label=\"Likelihood\" dp=\"4\" width=\"12\">"
	block.screen.log[length(block.screen.log)+1] <- "<likelihood idref=\"likelihood\"/>"
	block.screen.log[length(block.screen.log)+1] <- "</column>"
	block.screen.log[length(block.screen.log)+1] <- "<column label=\"rootHeight\" sf=\"6\" width=\"12\">"
	block.screen.log[length(block.screen.log)+1] <- "<parameter idref=\"treeModel.rootHeight\"/>"
	block.screen.log[length(block.screen.log)+1] <- "</column>"
	block.screen.log[length(block.screen.log)+1] <- "<column label=\"meanRate\" sf=\"6\" width=\"12\">"
	block.screen.log[length(block.screen.log)+1] <- "<parameter idref=\"meanRate\"/>"
	block.screen.log[length(block.screen.log)+1] <- "</column>"
	block.screen.log[length(block.screen.log)+1] <- "<column label=\"tmrca(root)\" sf=\"6\" width=\"12\">"
	block.screen.log[length(block.screen.log)+1] <- "<parameter idref=\"treeModel.rootHeight\"/>"
	block.screen.log[length(block.screen.log)+1] <- "</column>"
	if(mix.model==T){
		block.screen.log[length(block.screen.log)+1] <- "<column label=\"Index\" sf=\"6\" with=\"12\">"
		block.screen.log[length(block.screen.log)+1] <- "<rateStatistic idref=\"branchRates.distributionIndex\"/>"
		block.screen.log[length(block.screen.log)+1] <- "</column>"
	}
	block.screen.log[length(block.screen.log)+1] <- "</log>"
	
	block.res.log <- character()
	block.res.log[length(block.res.log)+1] <- paste("<log id=\"fileLog\" logEvery=\"" ,as.integer(mcmc.freq), "\" fileName=\"simulated.data.log\" overwrite=\"false\">", sep="")
	block.res.log[length(block.res.log)+1] <- "<posterior idref=\"posterior\"/>"
	block.res.log[length(block.res.log)+1] <- "<prior idref=\"prior\"/>"
	block.res.log[length(block.res.log)+1] <- "<likelihood idref=\"likelihood\"/>"
	block.res.log[length(block.res.log)+1] <- "<parameter idref=\"treeModel.rootHeight\"/>"	
	for(j in 1:n.cals){
		block.res.log[length(block.res.log)+1] <- paste("<tmrcaStatistic idref=\"tmrca(",cal.names[j], ")\"/>", sep="")
	}
	for(k in 1:2){
		block.res.log[length(block.res.log)+1] <- paste("<tmrcaStatistic idref=\"tmrca(", min.med.names[k], ")\"/>", sep="")
	}
		#Recording clock-specific statistics
	if(clock=="lognormal" && mix.model==F){
		block.res.log[length(block.res.log)+1] <- "<parameter idref=\"ucld.mean\"/>"
		block.res.log[length(block.res.log)+1] <- "<parameter idref=\"ucld.stdev\"/>"
		}else if(clock=="exponential" && mix.model==F){
				block.res.log[length(block.res.log)+1] <- "<parameter idref=\"uced.mean\"/>"
			}else if(mix.model==T){
				block.res.log[length(block.res.log)+1] <- "<parameter idref=\"ucld.mean\"/>"
				block.res.log[length(block.res.log)+1] <- "<parameter idref=\"ucld.stdev\"/>"
				block.res.log[length(block.res.log)+1] <- "<parameter idref=\"uced.mean\"/>"
				}
	block.res.log[length(block.res.log)+1] <- "<rateStatistic idref=\"meanRate\"/>"
	block.res.log[length(block.res.log)+1] <- "<rateStatistic idref=\"coefficientOfVariation\"/>"
	block.res.log[length(block.res.log)+1] <- "<rateCovarianceStatistic idref=\"covariance\"/>"
	if(mix.model==T){
		block.res.log[length(block.res.log)+1] <- "<parameter idref=\"branchRates.distributionIndex\"/>"
	}
	block.res.log[length(block.res.log)+1] <- "</log>"
	#Here would be the log trees block. This was not included for the simulations
	block.res.log[length(block.res.log)+1] <- "</mcmc>"
	if(mix.model==T){
		block.res.log[length(block.res.log)+1] <- "<mixtureModelLogAnalyser fileName=\"simulated.data.log\" burnin=\"200\" discreteVariable=\"branchRates.distributionIndex\"/>"
	}
	block.res.log[length(block.res.log)+1] <- "<report>"
	block.res.log[length(block.res.log)+1] <- "<property name=\"timer\">"
	block.res.log[length(block.res.log)+1] <- "<mcmc idref=\"mcmc\"/>"
	block.res.log[length(block.res.log)+1] <- "</property>"
	block.res.log[length(block.res.log)+1] <- "</report>"
	block.res.log[length(block.res.log)+1] <- "</beast>"

	all.blocks <- c(block.1, block.taxa, block.calib.tax, block.min.med, block.secs, block.pat.tree.model.size, block.tree.fix, block.tree.param, block.cal.stat, block.min.med.stat, block.spec, block.clock, block.model.ops, block.mcmc, block.screen.log, block.res.log)
	writeLines(all.blocks, "simulated.data.xml")
	#return(all.blocks)

}