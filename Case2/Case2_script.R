## ----------------------------------------------------------
## Script name: Case 1 example for simulating GEI
##
## Authors: Jon Bancic and Daniel Tolhurst
## Date Created: 2023-01-25
## Email: jbancic@ed.ac.uk
##
## ----------------------------------------------------------
## Description:
##  This script demonstrates the implementation of Case 1, 
##  when between-environment correlation is known.
## ----------------------------------------------------------
source(file = "../Simulation/Functions_asrGEI.R")

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

# OPTIONAL: heterogeneous genetic variance
varGe  = abs(diag(rnorm(n = p, mean = 1, sd = 0.5)))
Gmat   = sqrt(varGe) %*% Cmat %*% sqrt(varGe)

# 2-- Decompose G_e into p terms (without Cholesky)
#============================================================
Lam = svd(Gmat)$u
D   = diag(svd(Gmat)$d)

#-- Compare decomposed Gmat to simulated Cmat
range(Cmat - cov2cor(Lam %*% D %*% t(Lam)))
plot(cov2cor(Lam %*% D %*% t(Lam)), Cmat)
# Variance explained by terms
plot(diag(D)/sum(diag(D)), type = "line", 
     ylab = "Proportion of variance explained",
     xlab = "Environmental covariate",
     ylim = c(0,1))

# 2-- Optional: Decompose G_e into p terms (with Cholesky)
#============================================================
# Lam = svd(t(chol(Cmat)))$u
# D   = diag(svd(t(chol(Cmat)))$d^2)
# #-- Compare decomposed Gmat to simulated Cmat
# range(Cmat - cov2cor(Lam %*% D %*% t(Lam)))
# plot(cov2cor(Lam %*% D %*% t(Lam)), Cmat)
# # Variance explained by terms
# plot(diag(D)/sum(diag(D)), type = "line", 
#      ylab = "Proportion of variance explained",
#      xlab = "Environmental covariate",
#      ylim = c(0,1))

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

# 4-- OPTIONAL: reduced rank form with k terms
#============================================================
# Lam[,(k+1):p]    = 0
# diag(D)[(k+1):p] = 0

# 5-- Simulate genotype coefficients
#============================================================
# Scale and center coefficients to mean 0 and variance 1
f = scale(matrix(rnorm(v*p, mean = 0, sd = 1), ncol = p))
apply(f, 2, mean)
apply(f, 2, var)

# 6-- Construct the vector of genotypic values
#============================================================
varCoef = diag(rep(1, each = p)) # coefficient variance
GV      = matrix(f %*% t(sqrt(varCoef) %*% Lam %*% sqrt(D)), ncol = p)
#-- Compare Cmat based on GVs to simulated Cmat
range(cor(GV) - Cmat)
plot(cor(GV), Cmat)

# 7-- Construct the vector of phenotypic values
#============================================================
# Adjust number of environments, replcations and error variance 
# to obtain a desired heritability
VarA   = 1                 # additive genetic variance
mainG  = d[1]/sum(d)    
varA   = VarA*mainG        # variance due to main genetic effects  
varGEI = VarA*(VarA-mainG) # varinace due to non-crossover GEI    
varE   = 1
nEnvs  = 20
nReps  = 2
(h2    = varA/(varA + varGEI/nEnvs + varE/(nEnvs*nReps)))

#-- Simulate phenotypes
error = matrix(sqrt(varE/nReps)*rnorm(n = v*p, mean = 0, sd = 1), ncol = p)
sd(apply(error, 2, var))
pheno = GV + error

# OPTIONAL: heterozygous residual variance 
# varE <- diag(rnorm(p, mean = varE, sd = 0.5))
# error <- matrix(rnorm(n = v*p, mean = 0, sd = 1), ncol = p) %*% sqrt(varE/nReps)
# sd(apply(error, 2, var))
# pheno = GV + error

#-- Check heritability
# Mean-line heritability
cor(rowMeans(GV), rowMeans(pheno[,1:nEnvs]))^2
# Per environment heritability
diag(cor(GV[,1:nEnvs], pheno[,1:nEnvs]))^2       # mean per env
mean(diag(cor(GV[,1:nEnvs], pheno[,1:nEnvs])))^2 # per env
# Overall heritability
cor(c(pheno), c(GV))^2       
mean(diag(cor(pheno, GV))^2) # per env

# Compare true and phenotypic correlation matrices
plotCMAT(cmat = Cmat, order = T)
order <- plotCMAT(cmat = Cmat, order = T)$order
plotCMAT(cmat = cor(pheno), order = T, setorder = order)
range(cor(pheno) - Cmat)

#============================================================
## ----- Latent regression plots ----------------------------
temp <- data.frame(Env = rep(1:p,each=v),
                   aVar = as.factor(1:v),
                   Lam = (rep(Lam[,1],each=v)),
                   f = f[1:v])
# main effect
ggplot(temp, aes(y = f, x = Lam)) +
  geom_line(data=temp, aes(y = f, x = -Lam, group = aVar))

# main effect + environment effect
# ggplot(temp, aes(y = f, x = Lam)) +
#   geom_line(data=temp, aes(y = f, x = Lam, group = aVar))

# main effect + non-crossover
mean <- (-mean(temp$Lam))
ggplot(temp, aes(y = Lam*f, x = Lam)) +
  geom_vline(xintercept = mean, linetype = "dotted") +
  geom_line(data = temp, aes(y = Lam*f, x = -Lam, group = aVar)) +
  ggtitle("Factor 1")

# main effect + crossover
temp <- data.frame(Env  = rep(1:p,each=v),
                   aVar = as.factor(1:v),
                   Lam  = (rep(Lam[,1],each=v)),
                   Lam2 = (rep((Lam[,2]),each=v)),
                   Lam3 = (rep((Lam[,3]),each=v)),
                   f  = f[1:v],
                   f2 = f[201:400],
                   f3 = f[401:600])
ggplot(temp, aes(y = Lam2*f2, x = Lam)) +
  # geom_vline(xintercept = mean, linetype = "dotted") +
  geom_line(data=temp, aes(y = Lam2*f2, x = Lam2, group = aVar)) +
  ggtitle("Factor 2")
ggplot(temp, aes(y = Lam3*f3, x = Lam)) +
  # geom_vline(xintercept = mean, linetype = "dotted") +
  geom_line(data=temp, aes(y = Lam3*f3, x = Lam3, group = aVar)) +
  ggtitle("Factor 3")

## ----- OP-RMSD plot --------------------------------------------
# Overall performance is based on the first environmental covariate
OP <- mean(Lam[,1]) * f[1:v]    # Lam * scores
# RMSD = root mean squared deviations
dev  <- matrix(f, ncol = p)[,2:p] %*% t(Lam[,2:p]) # higher order factors
RMSD <- sqrt(rowMeans(dev^2))                      # square root of means
# Plot
df <- data.frame(OP,RMSD)
ggplot(df, aes(x=RMSD, y=OP)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0,linetype = "dotted") +
  geom_vline(xintercept = 0,linetype = "dotted") 
  
## ----- PCA biplot --------------------------------------------
colnames(pheno) <- paste("Env",1:dim(pheno)[2])
rownames(pheno) <- paste("Geno",1:dim(pheno)[1])
# perform PCA
results <- princomp(pheno)
# visualize results of PCA in biplot
biplot(results)