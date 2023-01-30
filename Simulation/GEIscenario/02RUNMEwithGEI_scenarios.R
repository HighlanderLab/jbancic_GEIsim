#!/exports/cmvm/eddie/eb/groups/HighlanderLab/communal/R-4.1.3/R-4.1.3/bin/Rscript
#$ -N GEIsim
#$ -cwd
#$ -R y
#$ -pe sharedmem 3
#$ -l h_vmem=8G
#$ -l h_rt=00:30:00
#$ -j y
#$ -V
#$ -P roslin_HighlanderLab
#$ -M jbancic@exseed.ed.ac.uk

#-- Load packages
rm(list = ls())
library(dplyr)
library(AlphaSimR)
source("../Functions_asrGEI.R")

#-- Load data
Rep = Sys.getenv("SGE_TASK_ID")
load(paste0("../../HAPLO", Rep, ".RData"))

#-- Import scenario-specific parameters
testScenario    = read.csv(file = "testScenario.csv",sep = " ",header = F)
weightMean      = as.integer(testScenario[1])/100
weightStability = as.integer(testScenario[2])/100

#-- Set scenario name
ScenarioName = paste0("withGEI_",testScenario[1],testScenario[2])

#-- Fill breeding pipeline with unique individuals from initial parents
source("FillPipeline_withGEI.R")

#-- Simulate a between environments and years correlation matrix 
source("CreateCorrelations.R")

#-- Initialize a placeholder for tracked variables
trackParams = data.frame(matrix(NA,nrow = 0,ncol = 16))

#-- Cycle years to make more advanced parents
j=0; n=0
for(year in 1:nCycles){ 
  time = timestamp(prefix = "",suffix = "",quiet = T)
  cat("Working on year: ",year,"    (",time,")  \n", sep = "")
  # Take correlations
  take    = ((nEnvs*(nCycles-year))+1):(nEnvs*(nCycles-year+1))
  tmpCmat = Cmat[take,take]
  source("../UpdateParents.R") 
  source("AdvanceYear_withGEI.R") 
}

#-- Assign names to variables
colnames(trackParams) = c("Year","Rep","Stage", "Scenario",
                          "MeanG", "VarG", "StabilityG",
                          "trueMean", "trueSD", "trueSI", 
                          "obsMean",  "obsSD",  "obsSI", 
                          "estMean",  "estSD",  "estSI")

#-- Save results
saveRDS(trackParams, paste0("trackParams_"  ,ScenarioName,".",Rep, ".rds"))
