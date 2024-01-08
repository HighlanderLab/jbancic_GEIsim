#!/exports/applications/apps/community/roslin/R/4.3.1-asreml/bin/Rscript
#$ -N GSGEI
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


# Script name: Genomic selection line breeding program with GEI
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# Description:
# This script considers a genomic selection line breeding program 
# with moderate GEI. Genomic model is used to predict GEBVs of HDRW 
# individuals and to improve selection in PYT and AYT stages.


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

#-- Selection criteria
weightMean = 1
weightSD   = 1-weightMean

#-- Set scenario name
Scenario = paste0("GSGEI_",weightMean,"-",weightSD)
subScenario  = "ModerateGEI"

#-- Load in pre-simulated TPE genetic variance matrix 
# (created in 00prepareTPE script)
df = readRDS(paste0("Ge_",subScenario,".rds"))

#-- Get environmental covariates by decomposing TPE genetic variance matrix 
aLamTPE = decompMat(df$aGe.TPE, kTerms)

#-- Obtain a set of 600 environments for simulation (same across reps)
aLamMET = aLamTPE[df$sample,]

# Reshuffle environments in aLamMET (different shuffle in each rep)
shuffle = sample(1:(nCycles*nEnvs), nCycles*nEnvs, replace = F)
aLamMET = aLamMET[shuffle,]
rownames(aLamMET) = colnames(df$aGe.TPE)[df$sample][shuffle]

# Order correlation matrix within each year and then reorder aLamMET accordingly
order = plotCmat(cor_mat = cov2cor(aLamMET%*%t(aLamMET)), den_order = T,
                 groups = as.list(as.data.frame(matrix(1:(nCycles*nEnvs),ncol = nCycles))))$order
aLamMET = aLamMET[order,]
# (1:1000)[df$sample][shuffle][order]

#-- Fill breeding pipeline with unique individuals from founder parents
source("FillPipeline_withGEI.R")

#-- Run 20 years of breeding
#=======================================================================
for(year in 1:nCycles){ 
  time = timestamp(prefix = "",suffix = "",quiet = T)
  cat("Working on year: ",year,"    (",time,")  \n", sep = "")
  
  #-- Obtain environmental covariates for that simulation year (or MET)
  if (year <= nYearsTP) {
    # Phenotypic selection first nYearsTP years where MET is just current year's data 
    METsample = c(matrix(((year-1)*nEnvs+1):((year+nYearsTP-1)*nEnvs),nrow=nEnvs))
    nEnvsMET = 20
  } else {
    # Reverse order so that current year's trials come first
    METsample = c(matrix(((year-nYearsTP)*nEnvs+1):(year*nEnvs),nrow=nEnvs)[,nYearsTP:1])
    nEnvsMET = nEnvs*nYearsTP # do not change
  }
  tmp_aLamMET = subsetLam(aLamMET, subset = METsample)
  
  if (year <= nYearsTP) {
    #-- Select new parents
    source("UpdateParents_withGEI.R")
    
    #-- Perform selections
    source("AdvanceYear_withGEI.R")
  } else {
    #-- Run models
    source("runGSModel.R") 
    
    #-- Select new parents
    source("UpdateParentsGS_withGEI.R")
    
    #-- Perform selections
    source("AdvanceYearGS_withGEI.R")
  }
  
  #-- Update training records
  source("StoreTrainPop.R")
}

#-- Condense results
source("condenseResults.R")

# Save results -----------------
cat("Saving rep:", Rep, "\n")
saveRDS(trackParams, paste0("Results_",subScenario,"_",Scenario,".",Rep,".rds"))
