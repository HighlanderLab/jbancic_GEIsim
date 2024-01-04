#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data
cat(" Advance year \n")

#-- Year 7
#Release variety
Release = selectIndGEI(pop = EYT, nInd = 1, use = selEYT)
# Report parameters
trackParams = rbind(trackParams,
                    reportParams(EYT,  year = year, Rep = Rep, 
                                 Pool = "NA",   Stage = "EYT", Scenario, subScenario,
                                 nEnvs = nEnvEYT, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET, aLamTPE = aLamTPE))

#-- Year 6
# Report parameters
trackParams = rbind(trackParams,
                    reportParams(AYT,  year = year, Rep = Rep, 
                                 Pool = "NA",   Stage = "AYT", Scenario, subScenario,
                                 nEnvs = nEnvAYT, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET, aLamTPE = aLamTPE))
# Perform selections - using genomic-assisted selection
EYT = selectIndGEI(pop = AYT, nInd = nEYT, use = selEYT)
EYT = setPhenoGEI(pop = EYT, aLam = tmp_aLamMET, 
                  nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT)
EYT = calcSelCriteria(pop = EYT, nEnvs = nEnvEYT,  nEnvsMET = nEnvsMET, aLamTPE = aLamTPE,
                      weightMean = weightMean, weightSD = weightSD)

#-- Year 5
# Report parameters
trackParams = rbind(trackParams,
                    reportParams(PYT,  year = year, Rep = Rep, 
                                 Pool = "NA",   Stage = "PYT", Scenario, subScenario,
                                 nEnvs = nEnvPYT, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET, aLamTPE = aLamTPE))
# Perform selections - using genomic-assisted selection
AYT = selectIndGEI(pop = PYT, nInd = nAYT, use = selAYT)
AYT = setPhenoGEI(pop = AYT, aLam = tmp_aLamMET,  
                  nEnvs = nEnvAYT, nReps = nRepAYT, varE = varE_AYT)
AYT = calcSelCriteria(pop = AYT, nEnvs = nEnvAYT,  nEnvsMET = nEnvsMET, aLamTPE = aLamTPE,
                      weightMean = weightMean, weightSD = weightSD)

#-- Year 4
# Report parameters
trackParams = rbind(trackParams,
                    reportParams(HDRW,  year = year, Rep = Rep, 
                                 Pool = "NA",   Stage = "HDRW", Scenario, subScenario,
                                 nEnvs = nEnvHDRW, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET, aLamTPE = aLamTPE))
# Perform selections - using genomic predicted values for selection
PYT = selectWithinFamGEI(HDRW, famMax, use = selPYT)
PYT = selectIndGEI(pop = PYT, nInd = nPYT, use = selPYT)

# Report parameters - check after selection in HDRW
trackParams = rbind(trackParams,
                    reportParams(PYT,  year = year, Rep = Rep, 
                                 Pool = "NA", Stage = "postHDRW", Scenario, subScenario,
                                 nEnvs = nEnvHDRW, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET, aLamTPE = aLamTPE))

PYT = setPhenoGEI(pop = PYT, aLam = tmp_aLamMET,  
                  nEnvs = nEnvPYT, nReps = nRepPYT, varE = varE_PYT)
PYT = calcSelCriteria(pop = PYT, nEnvs = nEnvPYT,  nEnvsMET = nEnvsMET, aLamTPE = aLamTPE,
                      weightMean = weightMean, weightSD = weightSD)

#-- Year 3
HDRW = setPhenoGEI(pop = DH, aLam = tmp_aLamMET,  
                   nEnvs = nEnvHDRW, nReps = nRepHDRW, varE = varE_HDRW)
HDRW = calcSelCriteria(pop = HDRW, nEnvs = nEnvHDRW,  nEnvsMET = nEnvsMET, aLamTPE = aLamTPE,
                       weightMean = weightMean, weightSD = weightSD)

#-- Year 2
DH = makeDH(F1, nDH)

#-- Year 1
F1 = randCross(Parents, nCrosses)