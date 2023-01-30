#===================================================================
#===================================================================
######### Modified ASR functions for simulations with GEI ##########
#===================================================================
#===================================================================

# Function sets phenotype based on the number environments and reps specified
#================================================================================
##' @param pop an object of Pop-class
##' @param Cmat between environment correlation matrix
##' @param nEnvs number of environments (p) in which genotypes are tested in. See details.
##' @param nReps number of replications within each environment. See details.
##' @param varE error (co)variances for traits. See details.
##' @param optionE specify 'hom' or 'het' to have homogeneous or heterogeneous error variance in different envrionments
##' @param k.terms specify the k (k < p) number of reduced rank terms to be used
##' @param mainG the proportion of variance explained by the main genetic effect (first term) and should range between 0 and 1
##' @param simParam an object of SimParam

setPhenoGEI = function(pop, Cmat = NULL, nEnvs = NULL, nReps = NULL, varE = NULL, 
                       optionE = "hom", k.terms = NULL, mainG = NULL, simParam = NULL) 
{
  match.arg(arg = optionE, 
            choices = c("noGEI","hom","het"))
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if (isSymmetric(Cmat) != TRUE) {
    stop("The correlation matrix is not symmetric")
  }
  if (nEnvs > 1 & is.null(optionE)) {
    stop("Specify if homogeneous ('hom') or heterogeneous ('het') error across environments")
  }
  if (nEnvs > dim(Cmat)[1]) {
    stop("nEnvs specified exceeds the dimensions of the provided Cmat")
  }
  if (nEnvs == 1) {optionE =  "singleEnv"}
  
  # Decompose correlation matrix
  nEnvTot = dim(Cmat)[1]
  Lam     = svd(t(chol(Cmat)))$u
  D       = diag(svd(t(chol(Cmat)))$d^2)
  if (!is.null(k.terms)) {
    if (nEnvs > dim(Cmat)[1]) {
      stop("k specified exceeds the dimensions of the provided Cmat")
    }
    if (!is.null(mainG)) {
      if (mainG < 0 | mainG > 1) {
        stop("mainG represents the proportion of variance explained by the main genetic effect,
           and should range between 0 and 1")
      }
      mainG      = mainG
      diag(D)[1] = mainG * nEnvTot
    } else {
      print("mainG is being calculated from the correlation matrix")
      d     = diag(D)
      mainG = d[1]/sum(d)
    }
    Lam[,k.terms:nEnvTot] = 0
    diag(D)[(k.terms+1):nEnvTot] = 0
    scale = (1-mainG)/sum(d[2:k.terms])*nEnvTot
    diag(D)[2:k.terms] = diag(D)[2:k.terms]*scale
    print(paste("Main G (first term) explains ",    round(diag(D)[1]/sum(D),2), 
                " and higher order terms explain ", round(sum(diag(D)[2:nTerms])/sum(D),2), 
                sep = ""))
  }
  varA    = diag(simParam$varA)
  GV      = matrix(pop@gv %*% t(sqrt(varA) %*% Lam %*% sqrt(D)), ncol = nEnvTot)

  if (optionE ==  "noGEI") {
    # When no GEI is considered
    if (any((Cmat[upper.tri(Cmat)] == 1) == FALSE)) {
      stop("Warning: When using 'noGEI' option, Cmat should have 1+1E-06 on diagonals and 1 on off-diagonals")
    }
    error  = matrix(rep(sqrt(varE/(nEnvs*nReps))*rnorm(n = pop@nInd, mean = 0, sd = 1), nEnvTot), ncol = nEnvTot)
  } else {
    # When single location and rep
    if (optionE ==  "singleEnv") {
      tmpGV = rowMeans(GV)
      error = matrix(sqrt(varE/nReps)*rnorm(n = pop@nInd, mean = 0, sd = 1), ncol = nEnvs)
    }
    # Simulate homogeneous error across locations 
    else if (optionE == "hom") {
      error = matrix(sqrt(varE/nReps)*rnorm(n = pop@nInd*nEnvTot, mean = 0, sd = 1), ncol = nEnvTot)
    }
    # Simulate heterogeneous error across locations
    else if (optionE == "het") {
      varE  = diag(rnorm(nEnvs, mean =  sqrt(varE), sd = 1))
      error = matrix(rnorm(n = pop@nInd*nEnvTot, mean = 0, sd = 1), ncol = nEnvTot) %*% sqrt(varE) 
    }
  }
  if (exists("tmpGV")) {
    pop@ebv       = GV
    pop@pheno[,1] = tmpGV + error
  } else {
    pop@ebv   = GV
    pop@pheno = GV + error
  }
  return(pop)
}
# Test
# pop = setPhenoGEI(pop = F1, Cmat = tmpCmat, nEnvs = 15, nReps = 1, varE = 4, optionE = "hom", simParam = SP)
# mean(diag(cor(pop@ebv, pop@pheno)^2))
# diag(cor(pop@ebv, pop@pheno)^2)

# Calculate OP and RMSD from phenotypic data
#=========================================================================
##' @param pop an object of Pop-class
##' @param nEnvs number of environments in which genotypes are tested in. See details.
##' @param weightMean selection index weight for mean genotype performance 
##' @param weightStability selection index weight for variability of a genotype

# when weightMean is set to 1, selection is based on phenotypic means only
# trueMean is a genotype genetic value calculated as an average of genetic values across all environments.
# obsMean is a genotype genetic value calculated as an average of genetic values across environments specific to a particular breeding stage.
# estMean is a genotype phenotypic value averaged across environments specific to a particular breeding stage.
# estSD is a genotype phenotypic standard deviation across environments specific to a particular breeding stage.
# estSI is a genotype selection index value given the weights assigned to OP and RMSD.
calcSelCriteria = function(pop, nEnvs = NULL, weightMean = 1, weightStability = 0) 
{
  if((weightMean+weightStability) != 1) {stop("Weights should sum to 1")}
  # True values
  trueMean = apply(pop@ebv, 1, mean)
  trueSD   = apply(pop@ebv, 1, sd)
  trueSI   = weightMean*trueMean - weightStability*trueSD
  pop = setMisc(pop, "trueMean", trueMean)
  pop = setMisc(pop, "trueSD",   trueSD)
  pop = setMisc(pop, "trueSI",   trueSI)
  # Observed values
  if (nEnvs == 1) {
    obsMean = pop@ebv[,1]
    obsSD   = NA
    obsSI   = NA
  } else {
    obsMean = apply(pop@ebv[,1:nEnvs], 1, mean)
    obsSD   = apply(pop@ebv[,1:nEnvs], 1, sd)
    obsSI   = weightMean*obsMean - weightStability*obsSD
  }
  pop = setMisc(pop, "obsMean", obsMean)
  pop = setMisc(pop, "obsSD",   obsSD)
  pop = setMisc(pop, "obsSI",   obsSI)
  # Estimated values
  if (nEnvs == 1) {
    estMean = pop@pheno[,1]
    estSD   = NA
    estSI   = NA
  } else {
    estMean = apply(pop@pheno[,1:nEnvs], 1, mean)
    estSD   = apply(pop@pheno[,1:nEnvs], 1, sd)
    estSI   = weightMean*estMean - weightStability*estSD
  }
  pop = setMisc(pop, "estMean", estMean)
  pop = setMisc(pop, "estSD",   estSD)
  pop = setMisc(pop, "estSI",   estSI)
  return(pop)
}
# Test
# pop = calcSelCriteria(pop = HDRW, nEnvs = nEnvHDRW,
#                 weightMean = weightMean, weightStability = weightStability)
# bind_rows(pop@misc)
# bind_rows(pop[1:5]@misc)
# bind_rows(pop[c(1,10)]@misc)

# Select individuals based on the desired measure 
#=================================================================
# use is any column in the "misc" slot to select genotypes
selectIndGEI = function(pop, nInd, use = "estMean", selectTop = TRUE)
{
  match.arg(arg = use, 
            choices = c("trueMean","trueSD","trueSI",
                        "obsMean","obsSD","obsSI",
                        "estMean","estSD","estSI"))
  response = matrix(unlist(getMisc(pop, use)), nrow = pop@nInd)
  if (is.matrix(response)) {
    stopifnot(ncol(response) == 1)
  }
  take = order(response, decreasing = selectTop)
  return(pop[take[1:nInd]])
}
# Test
# bind_rows(pop@misc)
# bind_rows(selectIndGEI(pop = pop, nInd = 5, use = "estMean")@misc)

# Select individuals based on the desired measure within family 
#=================================================================
selectWithinFamGEI = function(pop, nInd, use = "estMean", famType = "B", selectTop = TRUE)
{
  match.arg(arg = use, 
            choices = c("trueMean","trueSD","trueSI",
                        "obsMean","obsSD","obsSI",
                        "estMean","estSD","estSI"))
  families = AlphaSimR:::getFam(pop = pop, famType = "B")
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
#################### MET analysis functions ########################
#===================================================================
#===================================================================

# Function to analyse MET data and output BLUPs
#============================================================
##' @param pop is a correlation matrix
##' @param nEnvs is the variance
##' @param nReps

createMET = function(pop, nEnvs = NULL, nReps = NULL) 
{
  # Create MET dataset
  MET.data = data.frame(
    Genotype    = rep(pop@id, times = nEnvs),
    Environment = rep(paste("Env",1:nEnvs,sep = ""), each = nGenos),
    Rep         = rep(paste("Rep",1:nReps,sep = ""), each = nEnvs*nGenos),
    Phenotype   = c(pop@pheno))
  MET.data$Genotype = factor(MET.data$Genotype, 
                             levels = unique(as.character(MET.data$Genotype)))
  MET.data$Environment = factor(MET.data$Environment,
                                levels = unique(as.character(MET.data$Environment)))
  MET.data$Rep = factor(MET.data$Rep,
                        levels = unique(as.character(MET.data$Rep)))
  MET.data = MET.data[order(MET.data$Environment, MET.data$Rep, MET.data$Genotype),]
  # Create GRM
  M = pullSnpGeno(pop)-1
  K = mmat(x = M, center = F, scale = T, min.af = 0, max.mv = 1)
  rownames(K) = colnames(K) = MET.data$Genotype[1:nGenos]
  attr(K, "INVERSE") = FALSE
  diag(K) = diag(K) + 1e-6
  
  return(list(MET.data = MET.data, GRM = K))
}
# Test
# temp = createMET(pop = pop, nEnvs = 10, nReps = 2)
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

#===================================================================
#===================================================================
###################### Summary functions ###########################
#===================================================================
#===================================================================

# Extract values miscellaneous population parameters
#=========================================================================
getPopParams <- function(pop) 
{
  require(dplyr)
  return(bind_rows(pop@misc))
}
# Test
# mean(getPopParams(HDRW)$estMean)

# Mean genetic value and variance
#============================================================
summariseG = function(pop, Cmat) {
  Lam = svd(t(chol(Cmat)))$u
  D   = diag(svd(t(chol(Cmat)))$d^2)
  GV  = matrix(kronecker(rbind(colMeans(Lam)), diag(pop@nInd)) %*% c(pop@gv), ncol = 1)
  out = c(mean(GV), var(GV))
  names(out) = c("meanG", "varG")
  return(out)
}
# Test

# Calculate genetic mean, variance and variability of population
#=================================================================
trackGParams = function(pop) 
{
  # Gentic mean
  MeanG = mean(rowMeans(pop@ebv))
  # Genetic variance
  VarG  = mean(diag(var(pop@ebv)))
  # Genotype variability
  VariabG  = mean(apply(pop@ebv,1,sd))
  out = c(MeanG, VarG, VariabG)
  names(out) = c("MeanG", "VarG", "VariabG")
  return(out)
}
# Test
# trackGParams(pop = HDRW)

trackLambdas = function(pop, Cmat)
{
  Lam = svd(t(chol(Cmat)))$u
  D   = diag(svd(t(chol(Cmat)))$d^2)
  mean = colMeans(pop@gv %*% diag(colMeans(Lam %*% sqrt(D))))
  var  = apply(pop@gv %*% diag(colMeans(Lam %*% sqrt(D))),2,var)
  out = list(mean = mean, var = var)
  return(out)
}
# Test
# trackLambdas(EYT1, Cmat)
# pop@gv %*% t(sqrt(varA) %*% Lam %*% sqrt(D))
# colMeans(EYT1@gv %*% diag(colMeans(Lam %*% sqrt(D))))
# apply(EYT1@gv %*% diag(colMeans(Lam %*% sqrt(D))),2,var)


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


# Function for sampling environments from correlation matrix
#============================================================
# Either subset by sampling or just subseting first x
sampleCmat = function(Cmat, nEnvs = NULL, sample = FALSE)
{
  if (dim(Cmat)[1] < nEnvs) {
    stop("nEnvs is larger than correlation matrix dimension")
  }
  if (sample == TRUE) {
    take = sample(1:dim(Cmat)[1], nEnvs, replace = F)
    Cmat = Cmat[take,take]
  } else {
    Cmat = Cmat[1:nEnvs,1:nEnvs]
  }
  return(Cmat)
}
# Test
# sampleCmat(Cmat = Cmat[1:8,1:8], nEnvs = 5, sample = F)
# sampleCmat(Cmat = Cmat[1:8,1:8], nEnvs = 5, sample = T)


# Function plots and orders a correlation matrix
#===================================================================
plotCMAT <- function(cmat, order = F, setorder = NULL) 
{
  packages <- list("cluster","tidyr","ggplot2")
  do.call("require", packages)
  
  Cmata <- as.matrix(cmat)
  if (isSymmetric(Cmata) != TRUE) {stop("Matrix not symmetric")}
  if (order == T) {
    if (is.null(setorder)) {
      rownames(Cmata) <- colnames(Cmata) <- 1:dim(cmat)[1]
      Dmat <- 1 - Cmata
      agnes_e <- agnes(x = Dmat, diss = TRUE, method = "average")
      setorder <- agnes_e$order
      cor.df <- gather(as.data.frame(Cmata[setorder, setorder]))
      names(cor.df) <- c('Env1', "Cor")
      cor.df$Env1 <- factor(cor.df$Env1, levels = rev(setorder))
      cor.df$Env2 <- factor(rep(setorder, times = dim(cmat)[1]), levels = setorder)
      cor.df <- cor.df[order(cor.df$Env1, cor.df$Env2),]
      cor.df <- cor.df[,c("Env1","Env2","Cor")]
      cor.df$Cor[cor.df$Env1 == cor.df$Env2] <- NA
    } else {
      # stop(setorder != numeric)
      colnames(Cmata) <- rownames(Cmata) <- 1:dim(Cmata)[1]
      cor.df <- gather(as.data.frame(Cmata[setorder, setorder]))
      names(cor.df) <- c('Env1', "Cor")
      cor.df$Env1 <- factor(cor.df$Env1, levels = rev(setorder))
      cor.df$Env2 <- factor(rep(setorder, times = dim(cmat)[1]), levels = setorder)
      cor.df <- cor.df[order(cor.df$Env1, cor.df$Env2),]
      cor.df <- cor.df[,c("Env1","Env2","Cor")]
      cor.df$Cor[cor.df$Env1 == cor.df$Env2] <- NA
    }
  } 
  if (order == F) {
    colnames(Cmata) <- rownames(Cmata) <- 1:dim(Cmata)[1]
    cor.df <- gather(as.data.frame(Cmata))
    names(cor.df) <- c('Env1', "Cor")
    cor.df$Env1 <- factor(cor.df$Env1, levels = rev(colnames(Cmata)))
    cor.df$Env2 <- factor(rep(colnames(Cmata), times = dim(Cmata)[1]), levels = colnames(Cmata))
    cor.df <- cor.df[order(cor.df$Env1, cor.df$Env2),]
    cor.df <- cor.df[,c("Env1","Env2","Cor")]
    cor.df$Cor[cor.df$Env1 == cor.df$Env2] <- NA
  }
  hh <- rev(rainbow(256, start = 0, end = 2/3))
  # hh <- rev(rainbow(3))
  p <- ggplot(data = cor.df, aes(y= Env1, x = Env2, fill=Cor)) + geom_tile()+
    labs(x="Environment", y= "Environment") +
    scale_fill_gradientn(colours = hh, na.value = "white", limits= c(-1.1,1.1))
  
  temp <- list(plot = p)
  if (exists("setorder")) {
    temp <- append(temp,list(order = setorder))
  }
  return(temp)
}
# Test

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

