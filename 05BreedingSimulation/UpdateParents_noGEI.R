# Script removes 10 oldest parents each year and replaces them 
# with most recent genotypes from the EYT stage to form a new crossing block
cat(" Update parents \n")

# Drop 10 parents
Parents = Parents[11:nParents]
# Update with new 10 parents from the EYT stage
Parents = c(Parents, EYT)

#---------------------------------------------------------------------------------------------------------
#-- Assign new gvs and phenotypes to parents each year
Parents  = setPhenoGEI(pop = Parents, aLam = tmp_aLamMET,
                       nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT)
Parents  = calcSelCriteria(pop = Parents, nEnvs = nEnvEYT, nEnvsMET = nEnvsMET, 
                           weightMean = weightMean, weightSD = weightSD)
#---------------------------------------------------------------------------------------------------------

# Report parameters
trackParams = rbind(trackParams,
                    reportParams(Parents, year = year, Rep = Rep, 
                                 Pool = "NA", Stage = "Parents", Scenario, subScenario,
                                 nEnvs = nEnvEYT, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET))
