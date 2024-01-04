# Script name: Simulate TPE and generate samples
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# Description:
# Obtain genetic variance matrix for breeding simulation from 
# 04ModelComparison folder. Create a sample of 600 environments 
# which will be used over the course of a 20-year simulation period. 
# These 600 environments are shuffled in each simulation replicate.

# ---- Clean environment and load functions and packages ----

rm(list = ls())
source("Functions_asrGEI.R")

# Load TPE with moderate GEI from 04ModelComparison folder
load("../04ModelComparison/samplesModerateGEI.RData")

# Sample 600 environments 
sample <- sample(1:1000, 600, replace = F)
# **The same 600 envs are randomised each rep

# Take genetic variance matrix
aGe.TPE <- df$Ge
colnames(aGe.TPE) <- rownames(aGe.TPE) <- 1:dim(aGe.TPE)[1]

# Create data frame with Ge and samples 
df <- list("aGe.TPE" = aGe.TPE, 
           "sample"  = sample)

# Save data frame to be used in simulation with GEI
saveRDS(df, file = "Ge_ModerateGEI.rds")
