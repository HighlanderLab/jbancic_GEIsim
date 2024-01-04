#!/exports/applications/apps/community/roslin/R/4.3.1-asreml/bin/Rscript
#$ -N createHaplo
#$ -cwd
#$ -R y
#$ -pe sharedmem 1
#$ -l h_vmem=4G
#$ -l h_rt=00:10:00
#$ -j y
#$ -V
#$ -P roslin_HighlanderLab
#$ -M jbancic@exseed.ed.ac.uk

# Script name: Create founder haplotypes
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# Description:
# This script creates founder haplotypes for all scenarios.

# ---- Clean environment and load functions and packages ----

rm(list = ls())

# Load packages
packages <- c("AlphaSimR")
# install.packages(pkgs = packages, dependencies = TRUE)
lapply(packages, library, character.only = TRUE)

#-- Load data
Rep = Sys.getenv("SGE_TASK_ID")

#-- Load global parameters
source("GlobalParameters.R")

#-- Generate initial haplotypes
FOUNDERPOP = runMacs(nInd     = nParents, 
                     nChr     = 10, 
                     segSites = nQtl + nSnp,
                     inbred   = TRUE, 
                     species  = "WHEAT")

SP = SimParam$new(FOUNDERPOP)
SP$restrSegSites(nQtl,nSnp)

if(nSnp > 0){SP$addSnpChip(nSnp)}

# Simulate genotype slopes
SP$addTraitA(nQtlPerChr = nQtl,
             mean       = initMeanG,
             var        = initVarG, 
             corA       = corA)

Parents = newPop(FOUNDERPOP)
rm(FOUNDERPOP)

## Save R environment
save.image(paste0("HAPLO", Rep, ".RData"))
