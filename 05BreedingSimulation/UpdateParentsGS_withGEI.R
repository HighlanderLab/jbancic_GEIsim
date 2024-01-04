# Script removes 10 oldest parents each year and replaces them 
# with 10 genotypes from HDRW stage based on their predicted values
# to form a new crossing block
cat(" Update parents \n")

# Drop 10 parents
Parents = Parents[11:nParents]

# Update with new 10 parents from the EYT stage
newParents = selectWithinFamGEI(HDRW, nInd = 1, use = "estMean")
newParents = selectIndGEI(HDRW, nInd = 10, use = "estMean")

Parents = c(Parents, newParents)

#---------------------------------------------------------------------------------------------------------
#-- Assign new gvs and phenotypes to parents each year
Parents  = setPhenoGEI(pop = Parents, aLam = tmp_aLamMET,
                       nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT)
Parents  = calcSelCriteria(pop = Parents, nEnvs = nEnvEYT, nEnvsMET = nEnvsMET, aLamTPE = aLamTPE,
                           weightMean = weightMean, weightSD = weightSD)
#---------------------------------------------------------------------------------------------------------

# Report parameters
trackParams = rbind(trackParams,
                    reportParams(Parents, year = year, Rep = Rep, 
                                 Pool = "NA", Stage = "Parents", Scenario, subScenario,
                                 nEnvs = nEnvEYT, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET, aLamTPE = aLamTPE))
