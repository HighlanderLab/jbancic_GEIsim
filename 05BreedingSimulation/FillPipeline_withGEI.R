#-- Number of environments in MET
nEnvsMET = 20 # do not change

#-- Assign phenotypes with high accuracy to the initial parents
tmp_aLamMET = subsetLam(aLamMET, nEnvs = nEnvs*nYearsTP, year = 1)
Parents  = setPhenoGEI(pop = Parents, aLam = tmp_aLamMET,
                       nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT)
Parents  = calcSelCriteria(pop = Parents, nEnvs = nEnvEYT, nEnvsMET = nEnvsMET, aLamTPE = aLamTPE, 
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
    tmp_aLamMET = subsetLam(aLamMET, nEnvs = nEnvs*nYearsTP, 4)
    
    HDRW = setPhenoGEI(pop = DH, aLam = tmp_aLamMET, 
                       nEnvs = nEnvHDRW, nReps = nRepHDRW, varE = varE_HDRW)
    HDRW = calcSelCriteria(pop = HDRW, nEnvs = nEnvHDRW, nEnvsMET = nEnvsMET, aLamTPE = aLamTPE,  
                          weightMean = weightMean, weightSD = weightSD)
    cat(" HDRW h2 =", cor(HDRW@gv[,1], HDRW@pheno[,1])^2, "\n", sep = " ")
  }
  if(year<5){
    #Year 4
    tmp_aLamMET = subsetLam(aLamMET, nEnvs = nEnvs*nYearsTP, 3)
  
    PYT = selectWithinFamGEI(HDRW, famMax, use = "estMean")
    PYT = selectIndGEI(pop = PYT, nInd = nPYT, use = "estMean")
    PYT = setPhenoGEI(pop = PYT, aLam = tmp_aLamMET, 
                      nEnvs = nEnvPYT, nReps = nRepPYT, varE = varE_PYT)
    PYT = calcSelCriteria(pop = PYT, nEnvs = nEnvPYT, nEnvsMET = nEnvsMET, aLamTPE = aLamTPE, 
                          weightMean = weightMean, weightSD = weightSD)
    cat(" PYT h2  =", cor(rowMeans(PYT@gv[,1:nEnvPYT]), rowMeans(PYT@pheno[,1:nEnvPYT]))^2, "\n", sep = " ")
  }
  if(year<4){
    #Year 5
    tmp_aLamMET = subsetLam(aLamMET, nEnvs = nEnvs*nYearsTP, 2)
    
    AYT = selectIndGEI(pop = PYT, nInd = nAYT, use = "estMean")
    AYT = setPhenoGEI(pop = AYT, aLam = tmp_aLamMET, 
                      nEnvs = nEnvAYT, nReps = nRepAYT, varE = varE_AYT)
    AYT = calcSelCriteria(pop = AYT, nEnvs = nEnvAYT, nEnvsMET = nEnvsMET, aLamTPE = aLamTPE, 
                          weightMean = weightMean, weightSD = weightSD)
    cat(" AYT h2  =", cor(rowMeans(AYT@gv[,1:nEnvAYT]), rowMeans(AYT@pheno[,1:nEnvAYT]))^2, "\n", sep = " ")
  }
  if(year<3){
    #Year 6
    tmp_aLamMET = subsetLam(aLamMET, nEnvs = nEnvs*nYearsTP, 1)
    
    EYT = selectIndGEI(pop = AYT, nInd = nEYT, use = "estSI")
    EYT = setPhenoGEI(pop = EYT, aLam = tmp_aLamMET, 
                       nEnvs = nEnvEYT, nReps = nRepEYT, varE = varE_EYT)
    EYT = calcSelCriteria(pop = EYT, nEnvs = nEnvEYT, nEnvsMET = nEnvsMET, aLamTPE = aLamTPE,  
                           weightMean = weightMean, weightSD = weightSD)
    cat(" EYT h2  =", cor(rowMeans(EYT@gv[,1:nEnvEYT]), rowMeans(EYT@pheno[,1:nEnvEYT]))^2, "\n", sep = " ")
  }
  if(year<2){
    #Year 7
    Release = selectIndGEI(pop = EYT, nInd = 1, use = "estSI")
    cat(" Release", "\n")
  }
}
