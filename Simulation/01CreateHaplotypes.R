#!/exports/cmvm/eddie/eb/groups/HighlanderLab/communal/R-4.1.3/R-4.1.3/bin/Rscript
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

# Description: This script creates founder haplotypes for all scenarios

#-- Load packages
rm(list = ls())
library(AlphaSimR)

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

SP$addTraitAD(nQtlPerChr = nQtl,
              mean       = initMeanG,
              var        = initVarG, 
              meanDD     = initMeanD,
              varDD      = initVarD,
              corA       = corA,
              corD       = corD)

Parents = newPop(FOUNDERPOP)
rm(FOUNDERPOP)

#-- Add phenotype reflecting 2 years of evaluation in EYT
# Parents = setPheno(Parents, varE = varE_EYT, reps = nEnvEYT*nRepEYT, traits = 1)

## Save R environment
save.image(paste0("HAPLO", Rep, ".RData"))
