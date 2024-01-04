#===================================================================
#===================================================================
######### Modified ASR functions for simulations with GEI ##########
#===================================================================
#===================================================================

# Function sets phenotype based on the number environments and reps specified
#================================================================================
##' @param pop an object of Pop-class
##' @param aLam environmental covariates matrix with dimensions nEnvs x kTerms obtained from the decomposed between-environment additive correlation matrix
##' @param dLam environmental covariates matrix with dimensions nEnvs x kTerms obtained from the decomposed between-environment dominance correlation matrix
##' @param nEnvs number of environments (p) in which genotypes are tested 
##' @param nReps number of replications within each environment. If "noGEI" option is used, then nReps is the number of effective replications (nEnvs x nReps) 
##' and only the first environmental covariate is given a phenotype
##' @param errorType specify 'hom' or 'het' to have homogeneous or heterogeneous error variance across different envrionments. Alternatively,
##' specify 'noGEI' to obtain values based a single environment that can be adjusted by reps.
##' @param varE error (co)variances for environments. This value should decrease in the later stages of a breeding program where the size of plots and replications increases. 

setPhenoGEI = function(pop, aLam = NULL, dLam = NULL, nEnvs = NULL, nReps = NULL, varE = NULL, errorType = "hom") 
{
  match.arg(arg = errorType, 
            choices = c("noGEI","hom","het"))
  if (nEnvs > 1 & is.null(errorType)) {
    stop("Choose homogeneous ('hom') or heterogeneous ('het') error across environments")
  }
  tmp = setGenoGEI(pop, aLam, dLam, getScores = F)
  nEnvTot = dim(aLam)[1]
  # Add error
  if (errorType == "noGEI") {
    # When no GEI is considered - varE adjusted by effective replication
    error = matrix(rep(rnorm(pop@nInd, sd = sqrt(varE/(nEnvs*nReps))), nEnvTot), ncol = nEnvTot)
  } else if (errorType == "hom") {
    # Simulate homogeneous error across locations 
    error = matrix(rnorm(pop@nInd*nEnvTot, sd = sqrt(varE/nReps)), ncol = nEnvTot)
  } else if (errorType == "het") {
    # Simulate heterogeneous error across locations
    varE  = diag((rnorm(nEnvTot)*sqrt(0.01) + varE))
    error = matrix(rnorm(pop@nInd*nEnvTot), ncol = nEnvTot) %*% sqrt(varE/nReps)
  }
  pop@gv = tmp$gv
  pop@pheno = tmp$gv + error
  return(pop)
}
# Test
# tmp_dLam = NULL
# pop = setPhenoGEI(pop, aLam = tmp_aLamMET, dLam = tmp_dLam, nEnvs = 1, nReps = 1, varE = 8, errorType = "noGEI")
# mean(diag(cor(pop@gv, pop@pheno))^2)
# pop = setPhenoGEI(pop, aLam = tmp_aLamMET, dLam = tmp_dLam, nEnvs = 2, nReps = 1, varE = 4, errorType = "noGEI")
# mean(diag(cor(pop@gv, pop@pheno))^2)
# pop = setPhenoGEI(pop, aLam = tmp_aLamMET, dLam = tmp_dLam, nEnvs = 5, nReps = 2, varE = 4, errorType = "noGEI")
# mean(diag(cor(pop@gv, pop@pheno))^2)
# pop = setPhenoGEI(pop, aLam = tmp_aLamMET, dLam = tmp_dLam, nEnvs = 20, nReps = 2, varE = 4, errorType = "noGEI")
# mean(diag(cor(pop@gv, pop@pheno))^2)

# pop = setPhenoGEI(pop, aLam = tmp_aLamMET, nEnvs = 1, nReps = 1, varE = 12)
# cor(pop@gv[,1], pop@pheno[,1])^2
# pop = setPhenoGEI(pop, aLam = tmp_aLamMET, nEnvs = 2, nReps = 1, varE = 6)
# cor(rowMeans(pop@gv[,1:2]),rowMeans(pop@pheno[,1:2]))^2
# pop = setPhenoGEI(pop, aLam = tmp_aLamMET, nEnvs = 5, nReps = 2, varE = 6)
# cor(rowMeans(pop@gv[,1:5]),rowMeans(pop@pheno[,1:5]))^2
# pop = setPhenoGEI(pop, aLam = tmp_aLamMET, nEnvs = 20, nReps = 2, varE = 6)
# cor(rowMeans(pop@gv[,1:20]),rowMeans(pop@pheno[,1:20]))^2

# Function calculates additive, dominance and genetic values of a population
#================================================================================
setGenoGEI <- function(pop, aLam = NULL, dLam = NULL, simParam = NULL, getScores = FALSE) {
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  # if (exists("aLam")) {
  #   if (dim(aLam)[1] > pop@nTraits) {stop("aLam contains more environments than there are scores (i.e. nTraits)")} 
  # }
  nTraits <- pop@nTraits
  dP <- 1 # ploidy
  av <- dv <- gv <- NULL 
  QTL <- pullQtlGeno(pop, simParam = simParam) # nInd M nQTL
  for(i in seq_len(nTraits)){
    if (.hasSlot(simParam$traits[[i]], "domEff")) {
      av.i <- dv.i <-  gv.i <- NULL 
      # Recode QTL matrix
      A_m <- QTL - dP     # -1,0,1
      D_m <- QTL * (2-QTL) # 0,1,0
      # Genotypic Effects
      aEff <- simParam$traits[[i]]@addEff
      dEff <- simParam$traits[[i]]@domEff
      # Additive values
      av.i <- simParam$traits[[i]]@intercept + A_m %*% aEff
      # Dominance values
      dv.i <- D_m %*% dEff
      # Genetic values
      gv.i <- av.i + dv.i
      # Save values
      av <- cbind(av,rowSums(av.i))
      dv <- cbind(dv,rowSums(dv.i))
      gv <- cbind(gv,rowSums(gv.i))
    } else {
      av.i <- dv.i <-  gv.i <- NULL 
      # Recode QTL matrix
      A_m <- QTL - dP     # -1,0,1
      # Genotypic Effects
      aEff <- simParam$traits[[i]]@addEff
      # Additive values
      av.i <- simParam$traits[[i]]@intercept + A_m %*% aEff
      # Dominance values
      dv.i <- rep(0,dim(av.i)[1])
      # Genetic values
      gv.i <- av.i
      # Save values
      av <- cbind(av,rowSums(av.i))
      dv <- cbind(dv,dv.i)
      gv <- cbind(gv,rowSums(gv.i))
    }
  }
  if (getScores != TRUE) {
    nEnvTot = dim(aLam)[1]
    if (!is.null(dLam)) {
      # if (pop@nTraits != dim(dLam)[2]) {
      #   stop("dLam and scores (genotype coefficients ~ pop@nTraits) do not have the same number of terms")
      # }
      if (any((dim(aLam) == dim(dLam)) == FALSE)) {stop("Lam matrices are not of the same dimensions")}
      # dominance values 
      dv = matrix(dv %*% t(dLam), ncol = nEnvTot)
      # additive values
      av = matrix(av %*% t(aLam), ncol = nEnvTot)
      # total values
      gv = as.matrix(av + dv)
    } else {
      if (is.null(aLam)) {stop("Provide aLam if genetic values are desired")}
      gv = av = matrix(av %*% t(aLam), ncol = nEnvTot)
      dv = NA
    }
  }
  return(list("av" = av, "dv" = dv, "gv" = gv))
}
# Test
# ASR
# mean(apply(pop@gv,1,mean))
# var(apply(pop@gv,1,var))
# # my code
# mean(apply(setGenoGEI(pop, aLam = aLamMET)$gv,1,mean))
# var(apply(setGenoGEI(pop, aLam = aLamMET)$gv,1,var))

# Function decomposes correlation matrix to specified k number of terms
#================================================================================
##' @param mat provide between environment variance or correlation matrix
##' @param kTerms specify the k (k < p) number of reduced rank terms to be used
##' @param var specify variance within each environment

decompMat <- function(mat, kTerms = NULL, var = NULL) {
  if (isSymmetric(mat) != TRUE) {
    stop("Matrix is not symmetric")
  }
  S = -svd(mat)$u
  # if (as.logical(as.numeric(table(sign(S[,1]))[1]) > as.numeric(table(sign(S[,1]))[2]))) {S = -S}
  D = diag(sqrt(svd(mat)$d))
  if (!is.null(var)) {
    if (length(var) == 1) {
      var = diag(rep(var, dim(mat)[1]))
    }
    Lam = sqrt(var) %*% S %*% D
  } else {
    Lam = S %*% D
  }
  if (!is.null(kTerms)) {
    # reduced rank option
    if (kTerms > dim(mat)[1]) {
      stop("nTerms specified exceeds the dimensions of the provided matrix")
    }
    Lam = Lam[,1:kTerms]
  }
  return(Lam)
}
# Example
# mat <- df$aGe.TPE
# aLam <- decompMat(mat)
# range(df$aGe.TPE-aLam%*%t(aLam))
# 
# aLam <- decompMat(mat, kTerms = 60)
# range(df$aGe.TPE-aLam%*%t(aLam))
# 
# gv <- setGenoGEI(pop, aLam = aLam)$gv
# range(cov2cor(df$aGe.TPE)-cor(gv))
# mean(diag(cor(Cmat, cor(gv))))

# Function extracts lambdas given the year and the total number of environments
#================================================================================
subsetLam <- function(Lam, nEnvs, year, subset = NULL) {
  if (is.null(subset)) {
    if (year == 1) {
      Lam = Lam[1:nEnvs,]
    } else {
      Lam = Lam[((year-1)*nEnvs+1):(year*nEnvs),]
    }
  } else {
    if (length(subset) > dim(Lam)[1]) {stop("Subset larger than the row dimension of Lam")}
    Lam <- Lam[subset,]
  }
  return(Lam)
}
# Example
# tmp_aLam = subsetLam(Lam = aLamMET, nEnvs = nEnvs, year = 1)
# tmp_aLamMET = subsetLam(Lam = aLamTPE, subset = c(4,5,6,10,58))

# Calculate OP and RMSD from phenotypic data
#=========================================================================
##' @param pop an object of Pop-class
##' @param nEnvs number of environments in the current year's trials where genotypes are tested in. See details.
##' @param nEnvsMET number of environments in MET where genotypes were tested. See details.
##' @param aLamTPE decomposed TPE matrix obtained with decompMat() function to calculate true genetic values on the basis of a TPE. See details.
##' @param dLamTPE decomposed TPE matrix obtained with decompMat() function to calculate true genetic values on the basis of a TPE. See details.
##' @param weightMean selection index weight for mean genotype performance 
##' @param weightSD selection index weight for variability of a genotype

# when weightMean is set to 1, selection is based on phenotypic means only
# trueMean_TPE is a genotype genetic value calculated as an average of genetic values across all environments in a TPE
# trueMean_MET is a genotype genetic value calculated as an average of genetic values across all environments in a MET
# estMean is a genotype phenotypic value averaged across environments specific to a particular breeding stage in that year's MET
# estSD is a genotype phenotypic standard deviation across environments specific to a particular breeding stage in that year's MET
# estSI is a genotype selection index value given the weights assigned to OP and RMSD.
calcSelCriteria = function(pop, nEnvs = NULL, nEnvsMET = NULL, aLamTPE = NULL, dLamTPE = NULL, weightMean = 1, weightSD = 0) 
{
  if((weightMean+weightSD) != 1) {stop("Weights should sum to 1")}
  # True values based on TPE
  if (!is.null(aLamTPE)) {
    if (!is.null(dLamTPE)) {
      if (!identical(dim(aLamTPE), dim(dLamTPE))) {stop("aLam and dLam not of the same dimensions")}
    }
    tmp_pop = setPhenoGEI(pop, aLam = aLamTPE, dLam = dLamTPE, nEnvs = dim(aLamTPE)[1], nReps = 1, varE = 0)
    Mean = apply(tmp_pop@gv, 1, mean)
    tmp  = sweep(tmp_pop@gv, MARGIN = 2, STATS = colMeans(tmp_pop@gv)) # adjust by main environment effect
    SD   = apply(tmp, 1, sd)
    Mean_scaled = scale(apply(tmp, 1, mean), scale = T, center = F)
    SI  = weightMean*Mean_scaled - weightSD*SD
    pop = setMisc(pop, "Mean_TPE", Mean)
    pop = setMisc(pop, "SD_TPE",   SD)
    pop = setMisc(pop, "SI_TPE",   SI)
  } else {
    pop = setMisc(pop, "Mean_TPE", NA)
    pop = setMisc(pop, "SD_TPE",   NA)
    pop = setMisc(pop, "SI_TPE",   NA)
  }
  # True values based on MET
  if (nEnvsMET > dim(pop@gv)[2]) {stop(paste("nEnvsMET is larger than the total number of envs in a population"))}
  Mean = apply(pop@gv[,1:nEnvsMET], 1, mean)
  tmp  = sweep(pop@gv[,1:nEnvsMET], MARGIN = 2, STATS = colMeans(pop@gv[,1:nEnvsMET])) # adjust by main environment effect
  SD   = apply(tmp, 1, sd)
  Mean_scaled = scale(apply(tmp, 1, mean), scale = T, center = F)
  SI  = weightMean*Mean_scaled - weightSD*SD
  pop = setMisc(pop, "Mean_MET", Mean)
  pop = setMisc(pop, "SD_MET",   SD)
  pop = setMisc(pop, "SI_MET",   SI)
  # True values based current year's subset of trials
  if (nEnvs == 1) {
    Mean_subMET = pop@gv[,1]
    SD_subMET   = NA
    SI_subMET   = NA
  } else {
    Mean_subMET = apply(pop@gv[,1:nEnvs], 1, mean)
    SD_subMET   = apply(tmp[,1:nEnvs], 1, sd)
    Mean_subMET_scaled = scale(apply(tmp[,1:nEnvs], 1, mean), scale = T, center = F)
    SI_subMET   = weightMean*Mean_subMET_scaled - weightSD*SD_subMET
  }
  pop = setMisc(pop, "Mean_subMET", Mean_subMET)
  pop = setMisc(pop, "SD_subMET",   SD_subMET)
  pop = setMisc(pop, "SI_subMET",   SI_subMET)
  # Estimated values
  if (nEnvs == 1) {
    estMean = pop@pheno[,1]
    estSD   = NA
    estSI   = NA
  } else {
    estMean = apply(pop@pheno[,1:nEnvs], 1, mean)
    tmp = sweep(pop@pheno, MARGIN = 2, STATS = colMeans(pop@pheno))
    estSD   = apply(tmp[,1:nEnvs], 1, sd)
    estMean_scaled = scale(apply(tmp[,1:nEnvs], 1, mean), scale = T, center = F)
    estSI   = weightMean*estMean_scaled - weightSD*estSD
  }
  pop = setMisc(pop, "estMean", estMean)
  pop = setMisc(pop, "estSD",   estSD)
  pop = setMisc(pop, "estSI",   estSI)
  return(pop)
}
# Test
# pop <- calcSelCriteria(pop, nEnvs = 5, nEnvsMET = 60, aLamTPE = aLamTPE,
#                        weightMean = weightMean, weightSD = weightSD)
# getPopParams(pop)
# cor(getPopParams(pop)$Mean_TPE,getPopParams(pop)$Mean_MET)
# cor(getPopParams(pop)$Mean_TPE,getPopParams(pop)$estMean)
# cor(getPopParams(pop)$Mean_MET,getPopParams(pop)$Mean_subMET)


# Select individuals based on the desired measure 
#=================================================================
# use is any column in the "misc" slot to select genotypes
selectIndGEI = function(pop, nInd, use = "estMean", selectTop = TRUE)
{
  match.arg(arg = use, 
            choices = names(pop@misc[[1]]))
  response = matrix(unlist(getMisc(pop, use)), nrow = pop@nInd)
  if (is.matrix(response)) {
    stopifnot(ncol(response) == 1)
  }
  take = order(response, decreasing = selectTop)
  return(pop[take[1:nInd]])
}
# Test
# head(getPopParams(pop)$estMean[order(getPopParams(pop)$estMean, decreasing = T)])
# bind_rows(selectIndGEI(pop, nInd = 5, use = "estMean")@misc)

# Select individuals based on the desired measure within family 
#=================================================================
selectWithinFamGEI = function(pop, nInd, use = "estMean", famType = "B", selectTop = TRUE)
{
  match.arg(arg = use, 
            choices = names(pop@misc[[1]]))
  families = AlphaSimR:::getFam(pop = pop, famType = famType)
  response = matrix(unlist(getMisc(pop, use)), nrow = pop@nInd)
  if (is.matrix(response)) {
    stopifnot(ncol(response) == 1)
  }
  warn = FALSE
  selInFam = function(selFam) {
    index = which(families %in% selFam)
    y     = response[index]
    index = index[order(y, decreasing = selectTop)]
    # index = index[index %in% eligible]
    if (length(index) < nInd) {
      warn <<- TRUE
      return(index)
    }
    else {
      return(index[0:nInd])
    }
  }
  take = unlist(sapply(unique(families), selInFam))
  if (warn) {
    warning("One or more families are smaller than nInd")
  }
  return(pop[take])
}
# Test
# selectWithinFamGEI(pop = PYT,nInd = 1)
# AlphaSimR:::getFam(pop = selectWithinFamGEI(pop = PYT,nInd = 1), famType = "B")

#===================================================================
#===================================================================
###################### Summary functions ###########################
#===================================================================
#===================================================================

# Function to extract population parameters
#=========================================================================
getPopParams <- function(pop) {
  require(dplyr)
  return(bind_rows(pop@misc))
  # return(do.call("rbind", TrainPop@misc))
}
# Test
# getPopParams(pop)
# getPopParams(pop)$estMean

# Function to calculate additive, dominance and genetic mean and variance
#=========================================================================
##' @param pop an object of Pop-class
##' @param nEnvs number of environments in the current year's trials where genotypes are tested in. See details.
##' @param nEnvsMET number of environments in MET where genotypes were tested. See details.
##' @details MeanG/A/D and VarG/A/D are based on MET and subMeanG/A/D and subVarG/A/D is based on a subset of MET

calcPopParams = function(pop, nEnvs, nEnvsMET, aLam = NULL, dLam = NULL) {
  tmp = setGenoGEI(pop, aLam, dLam, getScores = F)
  if (nEnvs > nEnvsMET) {stop("nEnvs cannot be > nEnvsMET")}
  if (!is.null(dLam)) {
    if (pop@nInd == 1) { 
      # dominance values 
      dv    = tmp$dv[,1:nEnvsMET]
      MeanD = mean(dv)
      subMeanD = ifelse(nEnvs == 1, dv[1],  mean(dv[1:nEnvs]))
      VarD  = NA; subVarD = NA
      # additive values
      av    = tmp$av[,1:nEnvsMET]
      MeanA = mean(av)
      subMeanA = ifelse(nEnvs == 1, av[1],  mean(av[1:nEnvs]))
      VarA  = NA; subVarA = NA
      # total values
      gv    = tmp$gv[,1:nEnvsMET]
      MeanG = mean(gv)
      subMeanG = ifelse(nEnvs == 1, gv[1],  mean(gv[1:nEnvs]))
      VarMeanG  = NA; subVarMeanG = NA; VarG = NA; subVarG = NA
    } else {
      # dominance values 
      dv    = tmp$dv[,1:nEnvsMET]
      MeanD = mean(apply(dv, 1, mean))
      subMeanD = ifelse(nEnvs == 1, mean(dv[,1]),  mean(apply(dv[,1:nEnvs], 1, mean)))
      VarD  = mean(diag(var(dv)))
      subVarD  = ifelse(nEnvs == 1, var(dv[,1]), mean(diag(var(dv[,1:nEnvs]))))
      # additive values
      av    = tmp$av[,1:nEnvsMET]
      MeanA = mean(apply(av, 1, mean))
      subMeanA = ifelse(nEnvs == 1, mean(av[,1]),  mean(apply(av[,1:nEnvs], 1, mean)))
      VarA  = mean(diag(var(av)))
      subVarA  = ifelse(nEnvs == 1, var(av[,1]), mean(diag(var(av[,1:nEnvs]))))
      # total values
      gv    = tmp$gv[,1:nEnvsMET]
      MeanG = mean(apply(gv, 1, mean))
      subMeanG = ifelse(nEnvs == 1, mean(gv[,1]),  mean(apply(gv[,1:nEnvs], 1, mean)))
      VarMeanG  = var(apply(gv, 1, mean))
      subVarMeanG  = ifelse(nEnvs == 1, var(gv[,1]),  var(apply(gv[,1:nEnvs], 1, mean)))
      VarG  = mean(diag(var(gv)))
      subVarG = ifelse(nEnvs == 1, 
                       var(gv[,1]), mean(diag(var(gv[,1:nEnvs]))))   
    }
  } else {
    if (pop@nInd == 1) {
      MeanD = subMeanD = VarD = subVarD = NA
      gv    = tmp$gv[,1:nEnvsMET]
      MeanG = MeanA = mean(gv)
      subMeanG = subMeanA = ifelse(nEnvs == 1, gv[1],  mean(gv[1:nEnvs]))
      VarMeanG = NA; subVarMeanG = NA; VarG = VarA = NA; subVarG = subVarA = NA
    } else {
      MeanD = subMeanD = VarD = subVarD = NA
      gv    = tmp$gv[,1:nEnvsMET]
      MeanG = MeanA = mean(apply(gv, 1, mean))
      subMeanG = subMeanA = ifelse(nEnvs == 1, mean(gv[,1]), mean(apply(gv[,1:nEnvs], 1, mean)))
      VarMeanG  = var(apply(gv, 1, mean))
      subVarMeanG  = ifelse(nEnvs == 1, var(gv[,1]),  var(apply(gv[,1:nEnvs], 1, mean)))
      VarG  = VarA = mean(diag(var(gv)))
      subVarG = subVarA = ifelse(nEnvs == 1, 
                                 var(gv[,1]), mean(diag(var(gv[,1:nEnvs]))))
    }
  }
  output = c(MeanA, subMeanA, VarA, subVarA, 
             MeanD, subMeanD, VarD, subVarD, 
             MeanG, subMeanG, VarG, subVarG, 
             VarMeanG, subVarMeanG)
  names(output) = c("MeanA","MeanA_sub","VarA","VarA_sub",
                    "MeanD","MeanD_sub","VarD","VarD_sub", 
                    "MeanG","MeanG_sub","VarG","VarG_sub",
                    "VarMeanG", "VarMeanG_sub")
  return(output)
}
# Test
# calcPopParams(pop, nEnvs = dim(aLamTPE)[1], nEnvsMET = dim(aLamTPE)[1], aLam = aLamTPE, dLam = NULL)
# calcPopParams(pop, nEnvs = nEnvs, nEnvsMET = nEnvsMET, aLam = aLam, dLam = dLam)
# calcPopParams(pop, nEnvs = 5, nEnvsMET = nEnvsMET, aLam = tmp_aLamMET, dLam = NULL)

# Function to collect different progress parameters 
#================================================================================
reportParams = function(pop, year, Rep, Pool, Stage, Scenario, subScenario, 
                        nEnvs, nEnvsMET, aLam, dLam = NULL, aLamTPE = NULL, dLamTPE = NULL) {
  if (is.null(aLamTPE) & is.null(dLamTPE)) {
    output = data.frame(year  = year, 
                        Rep   = Rep, 
                        Pop   = Pool, 
                        Stage = Stage, 
                        Scenario = Scenario,
                        subScenario = subScenario,
                        matrix(NA,pop@nInd,ncol = 14,T),
                        matrix(calcPopParams(pop, nEnvs = nEnvs, nEnvsMET = nEnvsMET, 
                                             aLam = aLam, dLam = dLam),pop@nInd, ncol = 14, T),
                        ifelse(nEnvs == 1, var(pop@pheno[,1]), mean(diag(var(pop@pheno[,1:nEnvs])))),
                        bind_rows(pop@misc)) #<<
    colnames(output) = c("Year","Rep", "Pool","Stage", "Scenario", "subScenario",
                         "MeanA_TPE", "MeanA_subTPE", "VarA_TPE", "VarA_subTPE", 
                         "MeanD_TPE", "MeanD_subTPE", "VarD_TPE", "VarD_subTPE", 
                         "MeanG_TPE", "MeanG_subTPE", "VarG_TPE", "VarG_subTPE","VarMeanG_TPE", "VarMeanG_subTPE",
                         "MeanA_MET", "MeanA_subMET", "VarA_MET", "VarA_subMET", 
                         "MeanD_MET", "MeanD_subMET", "VarD_MET", "VarD_subMET", 
                         "MeanG_MET", "MeanG_subMET", "VarG_MET", "VarG_subMET","VarMeanG_MET", "VarMeanG_subMET", "estVar_MET",
                         "Mean_TPE",   "SD_TPE",   "SI_TPE",
                         "Mean_MET",   "SD_MET",   "SI_MET",
                         "Mean_subMET","SD_subMET","SI_subMET",
                         "estMean_MET","estSD_MET","estSI_MET")
  } else {
    output = data.frame(year  = year, 
                        Rep   = Rep, 
                        Pop   = Pool, 
                        Stage = Stage, 
                        Scenario = Scenario,
                        subScenario = subScenario,
                        matrix(calcPopParams(pop, nEnvs = dim(aLamTPE)[1], nEnvsMET = dim(aLamTPE)[1], 
                                             aLam = aLamTPE, dLam = dLamTPE),pop@nInd, ncol = 14, T),
                        matrix(calcPopParams(pop, nEnvs = nEnvs, nEnvsMET = nEnvsMET, 
                                             aLam = aLam, dLam = dLam),pop@nInd, ncol = 14, T),
                        ifelse(nEnvs == 1, var(pop@pheno[,1]), mean(diag(var(pop@pheno[,1:nEnvs])))),
                        bind_rows(pop@misc)) #<<
    colnames(output) = c("Year","Rep", "Pool","Stage", "Scenario", "subScenario",
                         "MeanA_TPE", "MeanA_subTPE", "VarA_TPE", "VarA_subTPE", 
                         "MeanD_TPE", "MeanD_subTPE", "VarD_TPE", "VarD_subTPE", 
                         "MeanG_TPE", "MeanG_subTPE", "VarG_TPE", "VarG_subTPE","VarMeanG_TPE", "VarMeanG_subTPE",
                         "MeanA_MET", "MeanA_subMET", "VarA_MET", "VarA_subMET", 
                         "MeanD_MET", "MeanD_subMET", "VarD_MET", "VarD_subMET", 
                         "MeanG_MET", "MeanG_subMET", "VarG_MET", "VarG_subMET","VarMeanG_MET", "VarMeanG_subMET", "estVar_MET",
                         "Mean_TPE",   "SD_TPE",   "SI_TPE",
                         "Mean_MET",   "SD_MET",   "SI_MET",
                         "Mean_subMET","SD_subMET","SI_subMET",
                         "estMean_MET","estSD_MET","estSI_MET")
  }
  return(output)
}
# Test
# reportParams(pop = Parents,  year = year, Rep = Rep, Pool = "NA", Stage = "Parents", Scenario, subScenario,
#              aLam = tmp_aLamMET, nEnvsMET = nEnvsMET, aLamTPE = aLamTPE, nEnvs = nEnvEYT)

# Function to calculate means and variances of genotype scores
#================================================================================
trackScores = function(pop, aLam = NULL, dLam = NULL)
{
  tmp = setGenoGEI(pop, aLam, dLam, getScores = T)
  output = list("av" = list(mean = apply(tmp$av,2,mean),
                            var  = apply(tmp$av,2,var)),
                "dv" = list(mean = apply(tmp$dv,2,mean), 
                            var  = apply(tmp$dv,2,var)),
                "gv" = list(mean = apply(tmp$gv,2,mean), 
                            var  = apply(tmp$gv,2,var)))
  return(output)
}
# Test
# trackScores(MaleInbredL1_TC)

#===================================================================
#===================================================================
###### Functions to simulate correlations and visualise them #######
#===================================================================
#===================================================================

# Function simulates correlations (adapted from Hardin et al. 2013)
#===================================================================
##' @param groups number of groups
##' @param size size of each group
##' @param rho baseline correlation in the kth group
##' @param delta baseline noise between group
##' @param epsilon maximum entry-wise random noise
##' @param eidim the dimension of the noise space
##' @param skew positive [-1,0] or negative [0,1] skew to the distribution of correlations

simCmat = function (groups, size, rho, delta, epsilon, eidim, skew)
{
  if (skew < -1 & skew > 1) {
    stop("The skew value should be between -1 and 1")
  }
  invert = FALSE
  if (skew < 0) { 
    invert = TRUE 
    skew = abs(skew)
  }
  
  ndim <- sum(size)
  bigcor <- matrix(rep(delta, ndim * ndim), ncol = ndim)
  for (i in 1:groups)  # fill up each block diagonal
  {
    cor <- matrix(rep(rho[i], size[i] * size[i]), ncol = size[i])
    if (i == 1) {bigcor[1:size[1], 1:size[1]] <- cor}
    if (i != 1) {bigcor[(sum(size[1:(i - 1)]) + 1):sum(size[1:i]), 
                        (sum(size[1:(i - 1)]) + 1):sum(size[1:i])] <- cor}
  }
  diag(bigcor) <- 1 - epsilon
  
  eivect <- c()
  for (i in 1:ndim) 
  {
    ei <- runif(eidim, -1, 1*skew)
    eivect <- cbind(eivect, sqrt(epsilon) * ei/sqrt(sum(ei^2)))
  }
  bigE <- t(eivect) %*% eivect
  cor.nz <- bigcor + bigE
  if (invert == TRUE) {
    cor.nz = 1 - cor.nz
  }
  return(cor.nz)
}
# Test
# Simulate a correlation matrix between environments
# nEnvPerYear <- 15
# rho         <- 0.5 # baseline correlation in each year
# epsilon     <- 0.99 - max(rho) # noise
# Cmat <- simCmat(groups  = 1, 
#                 size    = nEnvPerYear,
#                 rho     = rho,
#                 delta   = 0, 
#                 epsilon = epsilon,
#                 eidim   = 6, 
#                 skew    = 0.6)
# hist(Cmat[upper.tri(Cmat)])

# Function plots and orders a correlation matrix
#===================================================================
plotCmat <- function(cor_mat, 
                     den_order = TRUE, 
                     groups = NULL,
                     axis_title = "Environment"){
  
  require(cluster)
  require(tidyr)
  require(ggplot2)
  
  tot_envs <- ncol(cor_mat)
  env_order1 <- 1:tot_envs
  if(is.null(colnames(cor_mat)) | is.null(rownames(cor_mat))){colnames(cor_mat) <- rownames(cor_mat) <- env_order1}
  cor_mat <- as.matrix(cor_mat)
  if(isSymmetric(cor_mat) != TRUE){stop("Matrix is not symmetric")}
  Vmat = FALSE
  if (length(table(diag(cor_mat))) > 1) {
    Vmat = TRUE; print("Matrix does not have 1s on the diagonal. Assuming this is a variance matrix")
  }
  env_order2 <- colnames(cor_mat)
  if (!is.null(groups)){
    if(length(unlist(groups)) != tot_envs){stop("Groups specified do not conform with 'cor_mat'")}
    if(!is.list(groups)){stop("'groups' must be a list")}
    n_groups <- length(groups)
    n_envs <- unlist(lapply(groups, function(x) length(x)))
  }
  if(den_order){
    if (!is.null(groups)){
      dis_mat_list <- lapply(groups, function(x) 1 - cor_mat[x,x])
      env_order1 <- unname(unlist(lapply(dis_mat_list, function(x) colnames(x)[agnes(x = x, diss = TRUE, method = "average")$order])))
      env_order1 <- match(env_order1, env_order2)
      env_order2 <- colnames(cor_mat)[env_order1]
    }
    if(is.null(groups)){
      dis_mat <- 1 - cor_mat
      env_order1 <- agnes(x = dis_mat, diss = TRUE, method = "average")$order
      env_order2 <- colnames(cor_mat)[env_order1]
    } 
  }
  cor.df <- gather(as.data.frame(cor_mat[env_order2, env_order2]))
  names(cor.df) <- c('Env1', "Cor")
  cor.df$Env1 <- factor(cor.df$Env1, levels = rev(env_order2))
  cor.df$Env2 <- factor(rep(env_order2, times = tot_envs), levels = env_order2)
  cor.df <- cor.df[order(cor.df$Env1, cor.df$Env2),]
  rownames(cor.df) <- NULL
  cor.df <- cor.df[,c("Env1","Env2","Cor")]
  if (Vmat == FALSE) {cor.df$Cor[cor.df$Env1 == cor.df$Env2] <- NA}
  
  # hh <- rev(rainbow(256, start = 0, end = 2/3))
  hh <- rev(rainbow(256, start = 0, end = 2/3))[c(5:100,150:250)]
  pp <- ggplot(data = cor.df, aes(y= Env1, x = Env2, fill=Cor)) + geom_tile() +
    labs(x = axis_title, y = axis_title) +
    scale_fill_gradientn(colours = hh, na.value = "white", limits= c(-1.1,1.1)) +
    theme(axis.title.y = element_text(vjust = 1),
          axis.title.x = element_text(vjust = -1))
  if (Vmat == TRUE) {
    # hh <- hsv(seq(0.5,0.65,length.out = 100))
    hh <- hsv(seq(0.5,0.68,length.out = 200), 1, 0.9)
    pp <- pp + scale_fill_gradientn(colours = hh, na.value = "white") +
      labs(fill = "Var")
  }
  
  cor.df2 <- data.frame(Source = "All correlations",
                        Cor = cor.df[order(cor.df$Env1, rev(cor.df$Env2)),][lower.tri(cor_mat, diag= F),]$Cor)
  cor.df2$Mean <- mean(cor.df2$Cor)
  qq <- ggplot(data = cor.df2, aes(Cor)) + geom_histogram(bins = 20, color = "black", fill = "grey") +
    labs(x = paste0("Pair-wise correlations between ", tolower(axis_title),"s")) + 
    geom_vline(aes(xintercept = Mean), linetype = 2, colour = "steelblue", linewidth = 0.6) + 
    theme(axis.title.y = element_text(vjust = 2),
          axis.title.x = element_text(vjust = -1))
  
  if (!is.null(groups)){
    frames <- data.frame(Env1 = cumsum(c(0.5, n_envs))[1:n_groups],
                         Env2 =  cumsum(c(tot_envs + 0.5, -n_envs))[1:n_groups])
    pp <- pp + geom_rect(data = frames, linewidth = 0.5, fill = NA, colour = "black",
                         aes(xmin = Env1, xmax = Env1+n_envs,
                             ymin = Env2, ymax = Env2-n_envs))
    
    cor.df3 <- data.frame(Source = "Correlations within groups",
                          Cor = unlist(lapply(groups, function(x) cor_mat[x,x][upper.tri(cor_mat[x,x])])))
    cor.df3$Mean <- mean(cor.df3$Cor)
    cor_mat2 <- cor_mat
    for(i in 1:n_groups){cor_mat2[groups[[i]],groups[[i]]] <- NA}
    cor.df4 <- data.frame(Source = "Correlations between groups",
                          Cor = cor_mat2[upper.tri(cor_mat2)])
    cor.df4 <- cor.df4[!is.na(cor.df4$Cor),]
    cor.df4$Mean <- mean(cor.df4$Cor)
    cor.df5 <- rbind(cor.df2, cor.df3, cor.df4)
    cor.df5$Source <- factor(cor.df5$Source, levels = c("Correlations within groups", "Correlations between groups", "All correlations"))
    qq <- ggplot(data = cor.df5, aes(Cor)) + geom_histogram(bins = 20, color = "black", fill = "grey") +
      labs(x = paste0("Pair-wise correlations between ", tolower(axis_title),"s")) + 
      geom_vline(aes(xintercept = Mean), linetype = 2, colour = "steelblue", linewidth = 0.6) + 
      theme(axis.title.y = element_text(vjust = 2),
            axis.title.x = element_text(vjust = -1)) + facet_wrap(~Source)
  }
  temp <- list(heatmap = pp, hist = qq, data = cor.df, order = env_order1)
  return(temp)
}
# Test
# plotCmat(cor_mat = cov2cor(aGe), den_order = F)$heatmap # TPE
# plotCmat(cor_mat = cov2cor(aGe), den_order = T)$heatmap # TPE ordered
# plotCmat(cor_mat = cov2cor(aGe_sub), den_order = F)$heatmap # MET  
# plotCmat(cor_mat = cov2cor(aGe_sub),den_order = T,groups = as.list(as.data.frame(matrix(1:400,ncol = 20))))$heatmap # MET order within env.

# Function for sampling environments from correlation matrix
#============================================================
# Either subset by sampling or just subseting first x
sampleCmat = function(Cmat, nEnvs = NULL, sample = FALSE, replace = F)
{
  if (dim(Cmat)[1] < nEnvs) {
    stop("nEnvs is larger than correlation matrix dimension")
  }
  colnames(Cmat) <- rownames(Cmat) <- 1:dim(Cmat)[1]
  if (sample == TRUE) {
    take = sample(1:dim(Cmat)[1], nEnvs, replace = replace)
    Cmat = Cmat[take,take]
  } else {
    Cmat = Cmat[1:nEnvs,1:nEnvs]
  }
  return(Cmat)
}
# Test
# sampleCmat(Cmat = aGe[1:8,1:8], nEnvs = 5, sample = T)
# sampleCmat(Cmat = Cmat[1:8,1:8], nEnvs = 5, sample = T)


# Convert correlation matrix to covariance matrix
#===================================================================
##' @param C is a correlation matrix
##' @param S is the variance
cor2cov <- function(C,D){
  if (length(D) == 1) {
    D * C
  } else {
    diag(sqrt(D)) %*% C %*% diag(sqrt(D))
  }
}

# Reorder correlation matrix
#=======================================================================
##' @param cmat correlation matrix
##' @groups specify the number of groups

reorderCmat <- function(cmat, groups) {
  packages <- list("cluster")
  sapply(packages, require, character.only = TRUE)
  dim   <- dim(cmat)[1]
  size  <- dim/groups
  take  <- lapply(1:groups, function(x) matrix(1:dim, ncol = groups)[,x])
  order <- c()
  count <- 0
  for(i in 1:groups){
    dmat  <- 1 - cmat[take[[i]],take[[i]]]
    tmpOrder <- agnes(x = dmat, diss = TRUE, method = "average")$order
    tmpOrder <- tmpOrder + count
    order <- c(order, tmpOrder)
    count <- count + size
  }
  return(cmat[order,order])
}

#===================================================================
#===================================================================
#################### MET analysis functions ########################
#===================================================================
#===================================================================

# Function to store MET data for statistical analysis
#============================================================
##' @param pop is a correlation matrix
##' @param nEnvs is the variance
##' @param nReps
##' @param swapGender swap male parents for females and vice versa; useful when collecting hybrid records

createMET = function(pop, Stage, Year, nEnvs = NULL, nReps = NULL, swapGender = FALSE) 
{
  data     = data.frame(
    Year   = Year,
    Stage  = Stage,
    Loc    = rep(paste("Loc",1:nEnvs,sep = ""), each = pop@nInd),
    Rep    = nReps,
    Id     = rep(pop@id, times = nEnvs),
    Male   = rep(pop@father, times = nEnvs),
    Female = rep(pop@mother, times = nEnvs),
    Comb   = rep(paste0(pop@father,":",pop@mother), times = nEnvs),
    Pheno  = c(pop@pheno[,1:nEnvs]))
  if (swapGender == TRUE) {
    data$Male   = rep(pop@mother, times = nEnvs)
    data$Female = rep(pop@father, times = nEnvs)
  } 
  return(data)
}
# Test
# temp = createMET(pop = HybridL2_YT, Stage = "L2", Year = year, nEnvs = nEnvL2, nReps = nRepL2)
# str(temp)

# Analyse MET - compound symmetric model
#============================================================
analyseMET = function(MET.data, nEnvs = 1, nReps = 1, weightOP = 1, weightRMSD = 0)
{
  MET.data = temp$MET.data
  asreml.options(trace = FALSE)
  asr = asreml(Phenotype ~ Environment, 
               random   = ~ Genotype + idv(Environment):Genotype + 
                 diag(Environment):Rep,
               residual = ~ dsum(~ units | Environment),
               na.action = na.method("include"),
               data = MET.data)
  asr = update(asr)
  vars = summary(asr)$varcom
  
  # Heritability 
  VarG = vars[grep("Genotype", rownames(vars)),1][1]
  VarEnv = vars[grep("Genotype", rownames(vars)),1][2]
  VarE = vars[grep("!R", rownames(vars)),1]
  cat("Mean h2:", mean(round(VarG/(VarG+VarEnv+(VarE/nReps)),3)), sep = "\n")
  
  # Get genotype main effects and stability
  OP = asr$coefficients$random[grep("Environment",invert = T, rownames(asr$coefficients$random)),]
  BLUP = asr$coefficients$random[grep("Env.*Gen", rownames(asr$coefficients$random)),]
  RMSD = sqrt(apply(matrix(BLUP, ncol = nEnvs), 1, var))
  
  if((weightOP+weightRMSD) != 1) {stop("Weights should sum to 1")}
  output = (weightOP*OP + weightRMSD*RMSD)
  names(output) = levels(MET.data$Genotype)
  return(output)
}
# Test
# tmp = analyseMET(MET.data = temp, nEnvs = 10, nReps = 2)
# tmp = analyseMET(MET.data = temp, nEnvs = 10, nReps = 2,weightOP = 0.8, weightRMSD = 0.2)

################################################
measureGEI <- function(cov_mat,
                       prop = TRUE,
                       disentangle = FALSE,
                       groups = NULL,
                       old_measure = FALSE){
  if(!isSymmetric(cov_mat)) {stop("cov_mat is not symmetric")}
  if(!mean(cov_mat) > 0) {stop("cov_mat is not positive (semi) definite")}
  tot_var1 <- mean(diag(cov_mat))
  tot_var2 <- 1
  if(prop == F){
    tot_var1 <- 1
    tot_var2 <- mean(diag(cov_mat))}
  if (disentangle == FALSE) {
    Gvar   <- sum(cov_mat)/ncol(cov_mat)^2/tot_var1 
    # equivalent to: mean(cov_mat)/tot_var1
    GEIvar <- tot_var2 - Gvar
    var <- list(Gvar = Gvar, GEIvar = GEIvar)
  } 
  if (disentangle == T) {
    if (old_measure == F) {
      Gvar   <- sum(cov_mat)/ncol(cov_mat)^2/tot_var1
      Hetvar <- sum((sqrt(diag(cov_mat)) - mean(sqrt(diag(cov_mat))))^2)/ncol(cov_mat)/tot_var1 # het of variance
      nGEI <- Gvar + Hetvar
      cGEI <- tot_var2 - nGEI
    } else {
      # if(length(unique(sign(ss$u[,1]))) > 1){stop("First principal component captures cGEI, not just nGEI")}
      nGEI <- svd(cov_mat,nu = 1)$d[1]/ncol(cov_mat)/tot_var1
      cGEI <- tot_var2 - nGEI
    }
    var <- list(nGEI = nGEI, cGEI = cGEI)
  }
  return(var)
  # groups is a list with items in each group, can be names or indexes
  # initial check to see if all indexes are in cov_mat
  # if(!is.null(groups)){
  #   # within
  #
  #   # between
  #   bGvar <- sum(cov_mat)/ncol(cov_mat)/tot_var
  #   var <- list(#within = list(Gvar = wGvar, GEIvar = wGEIvar),
  #               #between = list(Gvar = bGvar, GEIvar = bGEIvar),
  #               overall = list(Gvar = Gvar, GEIvar = GEIvar))
  # }
}

