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

#-- Description
# This script does not consider GEI although but still uses GEI framework 
# Some alterations are required:
# 1. Set weightMean = 1 and weightStability = 0
# 2. Increase varE, so that it reflects the desired heritability
# 3. Set up a between-correlation matrix with 1+1E-06 on diagonals and 1s on off-diagonals

#-- Load packages
rm(list = ls())
library(dplyr)
library(AlphaSimR)
source("../Functions_asrGEI.R")

#-- Load data
Rep = Sys.getenv("SGE_TASK_ID")
load(paste0("../HAPLO", Rep, ".RData"))

#-- Import scenario-specific parameters
weightMean      = 1
weightStability = 0 # selection on stability not possible!
#Increase error variance
varE_HDRW = 14
varE_PYT = 8
varE_AYT = 5
varE_EYT = 5

#-- Set scenario name
ScenarioName = "noGEI"

#-- Fill breeding pipeline with unique individuals from initial parents
source("FillPipeline_noGEI.R")

#-- Simulate a between environments and years correlation matrix 
Cmat = matrix(1, nEnvs, nEnvs)
diag(Cmat) = diag(Cmat) + 1E-06
Cmat = kronecker(diag(nCycles),Cmat)

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
  source("AdvanceYear_noGEI.R") 
}

#-- Assign names to variables
colnames(trackParams) = c("Year","Rep","Stage", "Scenario",
                          "MeanG", "VarG", "StabilityG",
                          "trueMean", "trueSD", "trueSI", 
                          "obsMean",  "obsSD",  "obsSI", 
                          "estMean",  "estSD",  "estSI")
trackParams[,c("StabilityG","trueSD","trueSI","obsSD","obsSI","estSD","estSI")] <- NA

#-- Save results
saveRDS(trackParams, paste0("trackParams_", ScenarioName,".",Rep, ".rds"))
