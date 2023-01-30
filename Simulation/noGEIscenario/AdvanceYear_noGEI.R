#Advance breeding program by 1 year
#Works backwards through pipeline to avoid copying data

#-- Year 7
#Release variety
Release = selectIndGEI(pop = EYT, nInd = 1, use = "estSI")
# Report parameters
trackParams = rbind(trackParams,
                    data.frame(year, Rep, Stage = "Release", ScenarioName, 
                               matrix(trackGParams(Release),Release@nInd,3,T),
                               bind_rows(Release@misc)))

#-- Year 6
EYT = selectIndGEI(pop = AYT, nInd = nEYT, use = "estSI")
EYT = setPhenoGEI(pop = EYT, Cmat = tmpCmat,
                  nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT, optionE = "noGEI")
EYT = calcSelCriteria(pop = EYT, nEnvs = nEnvEYT,
                      weightMean = weightMean, weightStability = weightStability)
#Report parameters
trackParams = rbind(trackParams,
                    data.frame(year, Rep, Stage = "EYT", ScenarioName, 
                               matrix(trackGParams(EYT),EYT@nInd,3,T),
                               bind_rows(EYT@misc)))

#-- Year 5
AYT = selectIndGEI(pop = PYT, nInd = nAYT, use = "estMean")
AYT = setPhenoGEI(pop = AYT, Cmat = tmpCmat, 
                  nEnvs = nEnvAYT, nReps = nRepAYT, varE = varE_AYT, optionE = "noGEI")
AYT = calcSelCriteria(pop = AYT, nEnvs = nEnvAYT, 
                      weightMean = weightMean, weightStability = weightStability)
#Report parameters
trackParams = rbind(trackParams,
                    data.frame(year, Rep, Stage = "AYT", ScenarioName, 
                               matrix(trackGParams(AYT),AYT@nInd,3,T),
                               bind_rows(AYT@misc)))

#-- Year 4
PYT = selectWithinFamGEI(HDRW, famMax, use = "estMean")
PYT = selectIndGEI(pop = PYT, nInd = nPYT, use = "estMean")
PYT = setPhenoGEI(pop = PYT, Cmat = tmpCmat, 
                  nEnvs = nEnvPYT, nReps = nRepPYT, varE = varE_PYT, optionE = "noGEI")
PYT = calcSelCriteria(pop = PYT, nEnvs = nEnvPYT, 
                      weightMean = weightMean, weightStability = weightStability)
#Report parameters
trackParams = rbind(trackParams,
                    data.frame(year, Rep, Stage = "PYT", ScenarioName, 
                               matrix(trackGParams(PYT),PYT@nInd,3,T),
                               bind_rows(PYT@misc)))

#-- Year 3
HDRW = setPhenoGEI(pop = DH, Cmat = tmpCmat, 
                   nEnvs = nEnvHDRW, nReps = nRepHDRW, varE = varE_HDRW, optionE = "noGEI")
HDRW = calcSelCriteria(pop = HDRW, nEnvs = nRepHDRW, 
                       weightMean = weightMean, weightStability = weightStability)
#Report parameters
trackParams = rbind(trackParams,
                    data.frame(year, Rep, Stage = "HDRW", ScenarioName, 
                               matrix(trackGParams(HDRW),HDRW@nInd,3,T),
                               bind_rows(HDRW@misc)))

#-- Year 2
DH = makeDH(F1, nDH)

#-- Year 1
F1 = randCross(Parents, nCrosses)
