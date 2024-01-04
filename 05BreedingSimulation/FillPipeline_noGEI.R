#-- Assign phenotypes with high accuracy to the initial parents
Parents = setPhenoGEI(pop = Parents, aLam = tmp_aLamMET,
                      nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT, errorType = "noGEI")
Parents = calcSelCriteria(pop = Parents, nEnvs = nEnvEYT, nEnvsMET = nEnvsMET,
                          weightMean = weightMean, weightSD = weightSD)

#-- Set initial yield trials with unique individuals
for(year in 1:7){
  cat("FillPipeline year:", year, "of 7\n")
  if(year<8){
    #Year 1
    F1 = randCross(Parents, nCrosses)
  }
  if(year<7){
    #Year 2
    DH = makeDH(F1, nDH)
  }
  if(year<6){
    #Year 3
    HDRW = setPhenoGEI(pop = DH, aLam = tmp_aLamMET, 
                       nEnvs = nEnvHDRW, nReps = nRepHDRW, varE = varE_HDRW, errorType = "noGEI")
    HDRW = calcSelCriteria(pop = HDRW, nEnvs = nEnvHDRW, nEnvsMET = nEnvsMET,
                           weightMean = weightMean, weightSD = weightSD)
    cat(" HDRW h2 =", cor(HDRW@gv[,1], HDRW@pheno[,1])^2, "\n", sep = " ")
  }
  if(year<5){
    #Year 4
    PYT = selectWithinFamGEI(HDRW, famMax, use = "estMean")
    PYT = selectIndGEI(pop = PYT, nInd = nPYT, use = "estMean")
    PYT = setPhenoGEI(pop = PYT, aLam = tmp_aLamMET, 
                      nEnvs = nEnvPYT, nReps = nRepPYT, varE = varE_PYT, errorType = "noGEI")
    PYT = calcSelCriteria(pop = PYT, nEnvs = nEnvPYT, nEnvsMET = nEnvsMET,
                          weightMean = weightMean, weightSD = weightSD)
    cat(" PYT h2  =", cor(PYT@gv[,1], PYT@pheno[,1])^2, "\n", sep = " ")
  }
  if(year<4){
    #Year 5
    AYT = selectIndGEI(pop = PYT, nInd = nAYT, use = "estMean")
    AYT = setPhenoGEI(pop = AYT, aLam = tmp_aLamMET, 
                      nEnvs = nEnvAYT, nReps = nRepAYT, varE = varE_AYT, errorType = "noGEI")
    AYT = calcSelCriteria(pop = AYT, nEnvs = nEnvAYT, nEnvsMET = nEnvsMET,
                          weightMean = weightMean, weightSD = weightSD)
    cat(" AYT h2  =", cor(AYT@gv[,1], AYT@pheno[,1])^2, "\n", sep = " ")
  }
  if(year<3){
    #Year 6
    EYT = selectIndGEI(pop = AYT, nInd = nEYT, use = "estMean")
    EYT = setPhenoGEI(pop = EYT, aLam = tmp_aLamMET, 
                      nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT, errorType = "noGEI")
    EYT = calcSelCriteria(pop = EYT, nEnvs = nEnvEYT, nEnvsMET = nEnvsMET,
                          weightMean = weightMean, weightSD = weightSD)
    cat(" EYT h2  =", cor(EYT@gv[,1], EYT@pheno[,1])^2, "\n", sep = " ")
  }
  if(year<2){
    #Year 7
    Release = selectIndGEI(pop = EYT, nInd = 1, use = "estMean")
    cat(" Release", "\n")
  }
}
