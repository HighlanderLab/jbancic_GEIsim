#!/exports/applications/apps/community/roslin/R/4.3.1-asreml/bin/Rscript
#$ -N GSnoGEI
#$ -cwd
#$ -R y
#$ -pe sharedmem 1
#$ -l h_vmem=50G
#$ -l h_rt=48:00:00
#$ -j y
#$ -V
#$ -P roslin_HighlanderLab
#$ -M jbancic@exseed.ed.ac.uk
#$ -m eas


# Script name: Simulate genomic selection line breeding program without GEI
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# Description:
# This script does not consider GEI although but still uses GEI framework 
# Some alterations are required:
# 1. Set weightMean = 1 and weightSD = 0 since selection on stability is not possible
# 2. Tweak varE slightly, so that it is similar to heritability in GEI scenario
# 3. Set up a between-correlation matrix with 1s on diagonals and 0s on off-diagonals
# 4. Use setting errorType = "noGEI" when using function setPhenoGEI to indicate there is no GEI 


# ---- Clean environment and load functions and packages
rm(list = ls())
source("Functions_asrGEI.R")

# Load packages
packages <- c("AlphaSimR","tidyverse","asreml","ggpubr","dplyr",
              "ggplot2","reshape2","gridExtra","coda","ASRgenomics")
# install.packages(pkgs = packages, dependencies = TRUE)
lapply(packages, library, character.only = TRUE)

#-- Load data
Rep = Sys.getenv("SGE_TASK_ID")
load(paste0("HAPLO", Rep, ".RData"))

#-- Import scenario-specific parameters
# Selection criteria 
weightMean = 1 # do not change selection on stability not possible
weightSD   = 0 # do not change selection on stability not possible

#-- Set scenario name
Scenario = "GS_noGEI"
subScenario  = "ModerateGEI"

#-- Simulate a between-environment genetic variance matrix 
Cmat = matrix(1.474541, kTerms, kTerms) 
# NOTE: 1.474541 is the mean of simulated variance in GEI scenarios

#-- Decompose correlation matrix
tmp_aLamMET = decompMat(Cmat, kTerms)

#-- Number of environments in MET
nEnvsMET = 20 # do not change

#-- Fill breeding pipeline with unique individuals from initial parents
source("FillPipeline_noGEI.R")

#-- Run 20 years of breeding
#=======================================================================
for(year in 1:nCycles){ 
  time = timestamp(prefix = "",suffix = "",quiet = T)
  cat("Working on year: ",year,"    (",time,")  \n", sep = "")
  
  if (year <= nYearsTP) {
    #-- Number of environments in MET
    nEnvsMET = 20 # do not change
    
    #-- Select new parents
    source("UpdateParents_noGEI.R") 
    
    #-- Perform selections
    source("AdvanceYear_noGEI.R") 
  } else {
    #-- Number of environments in MET
    nEnvsMET = nEnvs*nYearsTP # do not change
    
    #-- Run models
    source("runGSModel_noGEI.R") 
    
    #-- Select new parents
    source("UpdateParentsGS_noGEI.R") 
    
    #-- Perform selections
    source("AdvanceYearGS_noGEI.R") 
  }
  
  #-- Update training records
  source("StoreTrainPop_noGEI.R")
}


#-- Condense results
source("condenseResults.R")

# Save results -----------------
cat("Saving rep:", Rep, "\n")
saveRDS(trackParams, paste0("Results_",subScenario,"_",Scenario,".",Rep,".rds"))
