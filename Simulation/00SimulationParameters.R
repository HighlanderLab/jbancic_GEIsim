## ----------------------------------------------------------
##  Theoretical - no GEI
## ----------------------------------------------------------
# HDRW h2 = 0.06666667
nEnvs  <- 1
nReps  <- 1
varA   <- 1
varE   <- 14
(h2 <- varA/(varA + varE/(nEnvs*nReps)))

# PYT h2 = 0.2
nEnvs  <- 2
nReps  <- 1
varA   <- 1
varE   <- 8
(h2 <- varA/(varA + varE/(nEnvs*nReps)))

# AYT h2 = 0.6666667
nEnvs  <- 5
nReps  <- 2
varA   <- 1
varE   <- 5
(h2 <- varA/(varA + varE/(nEnvs*nReps)))

# EYT2 h2 = 0.8888889
nEnvs  <- 20
nReps  <- 2
varA   <- 1
varE   <- 5
(h2 <- varA/(varA + varE/(nEnvs*nReps)))

## ----------------------------------------------------------
##  Theoretical - with GEI
## ----------------------------------------------------------
# Specify your starting parameters
VarA   <- 1
meanG  <- 0.6
varA   <- VarA*meanG       # 0.6 of VarA is due to main genetic effect
varGEI <- VarA*(1-meanG)   # 0.4 of VarA is due GEI

# HDRW h2 = 0.06666667
# Mean-line heritability
varE   <- 8
nEnvs  <- 1
nReps  <- 1
(h2 <- varA/(varA + varGEI/nEnvs + varE/(nEnvs*nReps)))
# Plot-level heritability
VarA/(VarA+varGEI+varE) 

# PYT h2 = 0.2142857
varE   <- 4
nEnvs  <- 2
nReps  <- 1
(h2 <- varA/(varA + varGEI/nEnvs + varE/(nEnvs*nReps)))

# AYT h2 = 0.6818182
varE   <- 2
nEnvs  <- 5
nReps  <- 2
(h2 <- varA/(varA + varGEI/nEnvs + varE/(nEnvs*nReps)))

# EYT2 h2 = 0.8955224
varE   <- 2
nEnvs  <- 20
nReps  <- 2
(h2 <- varA/(varA + varGEI/nEnvs + varE/(nEnvs*nReps)))
# Plot-level heritability
VarA/(VarA+varGEI+varE) 

## ----------------------------------------------------------
##  Simulated - with GEI
## ----------------------------------------------------------
source("Functions_asrGEI.R")
# Initial parameters
p  = 20  # number of environments
v  = 200 # number of genotypes
k  = 6   # number of environmental covariates

# 1-- Simulate between-environment correlation matrix
#============================================================
rho     = 0.5 # baseline correlation in each year
epsilon = 0.99 - max(rho) # noise
Cmat    = simCmat(groups  = 1, 
                  size    = p,
                  rho     = rho,
                  delta   = 0, 
                  epsilon = epsilon,
                  eidim   = 6, 
                  skew    = 0.6)
#-- Summary of Cmat
summary(Cmat[upper.tri(Cmat)])
hist(Cmat[upper.tri(Cmat)],
     main = "",
     ylab = "Density", 
     xlab = "Correlations")
plotCMAT(cmat = Cmat, order = T)

# 2-- Decompose G_e into p terms (without Cholesky)
#============================================================
Lam = svd(Gmat)$u
D   = diag(svd(Gmat)$d)

#-- Compare decomposed Cmat to the initial Cmat
cov2cor(Lam %*% D %*% t(Lam)) - Cmat
range(Cmat - cov2cor(Lam %*% D %*% t(Lam)))
# Variance explained by terms
plot(diag(D)/sum(diag(D)), type = "line", 
     ylab = "Proportion of variance explained",
     xlab = "Environmental covariate",
     ylim = c(0,1))

# 3-- Obtain genetic measures
#============================================================
# Total effects
(d = diag(D))
sum(d*colMeans(Lam)^2) * p        # variance

# Main genetic effects 
mean(Lam[,1])                     # term mean 
d[1]                              # variance
d[1]/sum(d)                       # proportion

# Non-crossover (scale) effects
nLam = Lam[,1] - mean(Lam[,1])
mean(nLam)                        # mean
d[1]*(t(nLam) %*% nLam)           # variance
d[1]*(t(nLam) %*% nLam)/sum(d)    # proportion

# Crossover effects
sum(d[2:p])                       # variance
sum(d[2:p])/sum(d)                # proportion

# 5-- Simulate genotype coefficients
#============================================================
# Scale and center coefficients to mean 0 and variance 1
f = scale(matrix(rnorm(v*p, mean = 0, sd = 1), ncol = p))

# 6-- Construct the vector of genotypic values
#============================================================
varCoef = diag(rep(1, each = p)) # coefficient variance
GV      = matrix(f %*% t(sqrt(varCoef) %*% Lam %*% sqrt(D)), ncol = p)
# Compare correlations between GVs in different envs to Cmat
range(cor(GV) - Cmat)
plot(cor(GV), Cmat)

# 7-- Construct the vector of phenotypic values
#============================================================
# Find varE for a desired heritability
VarA   = 1
mainG  = d[1]/sum(d)
varA   = VarA*mainG         # ~0.6 of VarA due to main genetic effect
varGEI = VarA*(VarA-mainG)  # ~0.4 of VarA due to GEI

#HDRW
varE   = 8
nEnvs  = 1
nReps  = 1
(h2    = varA/(varA + varGEI/nEnvs + varE/(nEnvs*nReps)))
error = matrix(sqrt(varE/nReps)*rnorm(n = v*p, mean = 0, sd = 1), ncol = p)
pheno = GV + error
# Mean-line heritability
cor(rowMeans(GV), pheno[,1])^2 # <--
# Per environment heritability
cor(GV[,1], pheno[,1])^2 # per env
# Overall heritability
cor(c(pheno), c(GV))^2 # <--       
mean(diag(cor(pheno, GV))^2) # per env # <--

#PYT
varE   = 4
nEnvs  = 2
nReps  = 1
(h2    = varA/(varA + varGEI/nEnvs + varE/(nEnvs*nReps)))
error = matrix(sqrt(varE/nReps)*rnorm(n = v*p, mean = 0, sd = 1), ncol = p)
pheno = GV + error
# Mean-line heritability
cor(rowMeans(GV), rowMeans(pheno[,1:nEnvs]))^2 # <--
# Per environment heritability
diag(cor(GV[,1:nEnvs], pheno[,1:nEnvs]))^2       # mean per env
mean(diag(cor(GV[,1:nEnvs], pheno[,1:nEnvs])))^2 # per env
# Overall heritability
cor(c(pheno), c(GV))^2       # <--
mean(diag(cor(pheno, GV))^2) # per env # <--

#AYT
varE   = 2
nEnvs  = 5
nReps  = 2
(h2    = varA/(varA + varGEI/nEnvs + varE/(nEnvs*nReps)))
error = matrix(sqrt(varE/nReps)*rnorm(n = v*p, mean = 0, sd = 1), ncol = p)
pheno = GV + error
# Mean-line heritability
cor(rowMeans(GV), rowMeans(pheno[,1:nEnvs]))^2 # <--
# Per environment heritability
diag(cor(GV[,1:nEnvs], pheno[,1:nEnvs]))^2       # mean per env
mean(diag(cor(GV[,1:nEnvs], pheno[,1:nEnvs])))^2 # per env
# Overall heritability
cor(c(pheno), c(GV))^2       
mean(diag(cor(pheno, GV))^2) # per env

#EYT
varE   = 2
nEnvs  = 20
nReps  = 2
(h2    = varA/(varA + varGEI/nEnvs + varE/(nEnvs*nReps)))
error = matrix(sqrt(varE/nReps)*rnorm(n = v*p, mean = 0, sd = 1), ncol = p)
pheno = GV + error
# Mean-line heritability
cor(rowMeans(GV), rowMeans(pheno[,1:nEnvs]))^2 # <--
# Per environment heritability
diag(cor(GV[,1:nEnvs], pheno[,1:nEnvs]))^2       # mean per env
mean(diag(cor(GV[,1:nEnvs], pheno[,1:nEnvs])))^2 # per env
# Overall heritability
cor(c(pheno), c(GV))^2       
mean(diag(cor(pheno, GV))^2) # per env
