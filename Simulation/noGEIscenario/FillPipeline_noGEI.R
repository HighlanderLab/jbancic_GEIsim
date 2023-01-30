#-- Simulate a correlation matrix between environments
Cmat       = matrix(1, nEnvs, nEnvs)
diag(Cmat) = diag(Cmat) + 1E-06
Cmat       = kronecker(diag(6),Cmat)

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
    take = ((nEnvs*(5-year))+1):(nEnvs*(6-year))
    tmpCmat = Cmat[take,take]
    
    HDRW = setPhenoGEI(pop = DH, Cmat = tmpCmat, 
                       nEnvs = nEnvHDRW, nReps = nRepHDRW, varE = varE_HDRW, optionE = "noGEI")
    HDRW = calcSelCriteria(pop = HDRW, nEnvs = nRepHDRW, 
                           weightMean = weightMean, weightStability = weightStability)
    cat(" HDRW h2 =", cor(getPopParams(HDRW)$trueMean,getPopParams(HDRW)$estMean)^2, "\n", sep = " ")
  }
  if(year<5){
    #Year 4
    take = ((nEnvs*(4-year))+1):(nEnvs*(5-year))
    tmpCmat = Cmat[take,take]
    
    PYT = selectWithinFamGEI(HDRW, famMax, use = "estMean")
    PYT = selectIndGEI(pop = PYT, nInd = nPYT, use = "estMean")
    PYT = setPhenoGEI(pop = PYT, Cmat = tmpCmat, 
                      nEnvs = nEnvPYT, nReps = nRepPYT, varE = varE_PYT, optionE = "noGEI")
    PYT = calcSelCriteria(pop = PYT, nEnvs = nEnvPYT, 
                          weightMean = weightMean, weightStability = weightStability)
    cat(" PYT  h2 =", cor(getPopParams(PYT)$trueMean,getPopParams(PYT)$estMean)^2, "\n", sep = " ")
  }
  if(year<4){
    #Year 5
    take = ((nEnvs*(3-year))+1):(nEnvs*(4-year))
    tmpCmat = Cmat[take,take]
    
    AYT = selectIndGEI(pop = PYT, nInd = nAYT, use = "estMean")
    AYT = setPhenoGEI(pop = AYT, Cmat = tmpCmat, 
                      nEnvs = nEnvAYT, nReps = nRepAYT, varE = varE_AYT, optionE = "noGEI")
    AYT = calcSelCriteria(pop = AYT, nEnvs = nEnvAYT, 
                          weightMean = weightMean, weightStability = weightStability)
    cat(" AYT  h2 =", cor(getPopParams(AYT)$trueMean,getPopParams(AYT)$estMean)^2, "\n", sep = " ")
  }
  if(year<3){
    #Year 6
    take = ((nEnvs*(2-year))+1):(nEnvs*(3-year))
    tmpCmat = Cmat[take,take]
    
    EYT = selectIndGEI(pop = AYT, nInd = nEYT, use = "estSI")
    EYT = setPhenoGEI(pop = EYT, Cmat = tmpCmat, 
                      nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT, optionE = "noGEI")
    EYT = calcSelCriteria(pop = EYT, nEnvs = nEnvEYT, 
                          weightMean = weightMean, weightStability = weightStability)
    cat(" EYT  h2 =", cor(getPopParams(EYT)$trueMean,getPopParams(EYT)$estMean)^2, "\n", sep = " ")
  }
  if(year<2){
    #Year 7
    Release = selectIndGEI(pop = EYT, nInd = 1, use = "estSI")
    cat(" Release", "\n")
  }
}
