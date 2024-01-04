#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data
cat(" Advance year \n")

#-- Year 7
#Release variety
Release = selectIndGEI(pop = EYT, nInd = 1, use = "estMean")
# Report parameters
trackParams = rbind(trackParams,
                    reportParams(EYT,  year = year, Rep = Rep, 
                                 Pool = "NA",   Stage = "EYT", Scenario, subScenario,
                                 nEnvs = nEnvEYT, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET),
                    reportParams(Release,  year = year, Rep = Rep, 
                                 Pool = "NA",   Stage = "Release", Scenario, subScenario,
                                 nEnvs = nEnvEYT, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET))

#-- Year 6
# Report parameters
trackParams = rbind(trackParams,
                    reportParams(AYT,  year = year, Rep = Rep, 
                                 Pool = "NA",   Stage = "AYT", Scenario, subScenario,
                                 nEnvs = nEnvAYT, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET))
# Perform selections
EYT = selectIndGEI(pop = AYT, nInd = nEYT, use = "estMean")
EYT = setPhenoGEI(pop = EYT, aLam = tmp_aLamMET, 
                  nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT, errorType = "noGEI")
EYT = calcSelCriteria(pop = EYT, nEnvs = nEnvEYT, nEnvsMET = nEnvsMET,
                      weightMean = weightMean, weightSD = weightSD)

#-- Year 5
# Report parameters
trackParams = rbind(trackParams,
                    reportParams(PYT,  year = year, Rep = Rep, 
                                 Pool = "NA",   Stage = "PYT", Scenario, subScenario,
                                 nEnvs = nEnvPYT, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET))
# Perform selections
AYT = selectIndGEI(pop = PYT, nInd = nAYT, use = "estMean")
AYT = setPhenoGEI(pop = AYT, aLam = tmp_aLamMET,
                  nEnvs = nEnvAYT, nReps = nRepAYT, varE = varE_AYT, errorType = "noGEI")
AYT = calcSelCriteria(pop = AYT, nEnvs = nEnvAYT, nEnvsMET = nEnvsMET,
                      weightMean = weightMean, weightSD = weightSD)

#-- Year 4
# Report parameters
trackParams = rbind(trackParams,
                    reportParams(HDRW,  year = year, Rep = Rep, 
                                 Pool = "NA",   Stage = "HDRW", Scenario, subScenario,
                                 nEnvs = nEnvHDRW, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET))
# Perform selections
PYT = selectWithinFamGEI(HDRW, famMax, use = "estMean")
PYT = selectIndGEI(pop = PYT, nInd = nPYT, use = "estMean")
PYT = setPhenoGEI(pop = PYT, aLam = tmp_aLamMET,  
                  nEnvs = nEnvPYT, nReps = nRepPYT, varE = varE_PYT, errorType = "noGEI")
PYT = calcSelCriteria(pop = PYT, nEnvs = nEnvPYT, nEnvsMET = nEnvsMET,
                      weightMean = weightMean, weightSD = weightSD)

#-- Year 3
HDRW = setPhenoGEI(pop = DH, aLam = tmp_aLamMET,  
                   nEnvs = nEnvHDRW, nReps = nRepHDRW, varE = varE_HDRW, errorType = "noGEI")
HDRW = calcSelCriteria(pop = HDRW, nEnvs = nEnvHDRW, nEnvsMET = nEnvsMET,
                       weightMean = weightMean, weightSD = weightSD)

#-- Year 2
DH = makeDH(F1, nDH)

#-- Year 1
F1 = randCross(Parents, nCrosses)
