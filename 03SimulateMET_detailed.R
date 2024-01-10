# Script name: A framework to simulate phenotypic data based on a LMM
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# Description:
# This script demonstrates simulation of phenotypic data to construct
# a MET dataset. Two approaches are used to simulate the genetic effects:
#
# Approach 1: independent genotype slopes,
# Approach 2: correlated genotype slopes with additive effects simulated 
# using the AlphaSimR package.
#
# In both methods, plot error with spatial variation is simulated using 
# using the FieldSimR package and dataset is randomised using the
# agricolae package.


# ---- Clean environment and load functions and packages ----

rm(list = ls())
source("00SimulateFunctions.R")

# Load packages
packages <- c("AlphaSimR", "FieldSimR", "agricolae", 
              "asreml", "ASRgenomics")
# install.packages(pkgs = packages, dependencies = TRUE)
lapply(packages, library, character.only = TRUE)


# ---- Specify parameters ----

set.seed(123)
p  = 10   # No. environments
b  = 2    # No. replicated trials per envrionment 
v  = 200  # No. genotypes
mu = 4    # Overall trait mean
H2 = 0.3  # Plot-level heritability
k  = 7    # No. of multiplicative terms


########################################################################
# Approach 1: Simulate genotype slopes outside of AlphaSimR
########################################################################

# --- Simulate environmental effects ---
X   = kronecker(diag(p), rep(1, v * b))
tau = scale(rnorm(p))

# --- Simulate genetic effects ---
# 1. Simulate between-environment correlation matrix
Ce = simCmat(n_envs = p, mean_cor = 0.999 - 0.5, rank = 6, epsilon = 0.5)
plotCmat(Ce, den_order = TRUE)$heat
plotCmat(Ce, den_order = TRUE)$hist

# 2. Obtain and decompose between-environment variance matrix
# and take k terms
De = diag(1/rgamma(p, shape = 5, rate = 5))
Ge = sqrt(De) %*% Ce %*% sqrt(De)
U = svd(Ge)$u[,1:k]
L = diag(svd(Ge)$d[1:k])

# 3. Obtain covariates and slopes
S = U
slopes = scale(matrix(rnorm(k * v), ncol = k))
slopes = c(matrix(slopes, ncol = k) %*% sqrt(L))

# 4. Construct genotype by environment interaction effects
u = kronecker(S, diag(v)) %*% slopes

# Obtain genotype main effects
ug = rowMeans(matrix(u, ncol = p))

# Construct an initial MET data frame
df.MET = data.frame(id  = factor(rep(1:v, each = b)),
                    env = factor(rep(1:p, each = v * b)),
                    rep = factor(1:b),
                    u   = rep(u, each = b))
df.MET = df.MET[order(df.MET$env, df.MET$rep),]

# --- Assign Randomized Complete Block Design with agricolae ---
design = design.rcbd(trt = 1:v, r = p * b, seed = 0)$book

# Reorder by trial
take = matrix(1:(v * b * p), ncol = v, byrow = TRUE)
for(i in 1:(p * b)){
  df.MET[take[i, ], ] = df.MET[take[i,],][design$`1:v`[take[i,]],]
}

# TO DO: add unbalanced design

# --- Simulate plot errors and spatial variation with FieldSimR ---
H2    = abs(rnorm(p, H2, 0.1))
H2[H2 < 0] = 0; H2[H2 > 1] = 1
var_R = diag(De) / H2 - diag(De)
df.error = field_trial_error(n_envs = p,
                             n_blocks = b,
                             n_traits = 1,
                             n_cols = 20,
                             n_rows = 20,
                             var_R  = var_R,
                             prop_spatial = 0.1)

# Plot errors
e = df.error$e.Trait.1

# --- Create phenotypes ---
y = mu + X %*% tau + df.MET$u + e

# --- Construct MET dataset ---
df.MET = data.frame(env = df.MET$env,
                    rep = df.MET$rep,
                    col = df.error$col,
                    row = df.error$row,
                    id  = df.MET$id,
                    y   = y,
                    u   = df.MET$u,
                    e   = e)
head(df.MET)

# --- Run model ---
asr = asreml(y ~ 1 + env,
             random    = ~ rr(env, 4):id + diag(env):id,
             residual  = ~ dsum(~ar1(col):ar1(row) | env),
             # residual  = ~ dsum(~units|env),
             na.action = na.method("include"),
             data = df.MET)

while (!asr$converge) {
  asr = update.asreml(asr)
}

# Report
summary(asr)$varcomp
estInt.reg = matrix(asr$coefficients$random[grep("^rr.*id", 
        rownames(asr$coefficients$random))][1:(v * p)],ncol = p,byrow = F)
estInt.diag = matrix(asr$coefficients$random[grep("^env.*id", 
        rownames(asr$coefficients$random))][1:(v * p)],ncol = p,byrow = F)
Lam = matrix(asr$vparameters[grep("rr.*fa", names(asr$vparameters))], ncol = 4)
Psi = diag(asr$vparameters[grep("^env.*id", names(asr$vparameters))])
estR = asr$vparameters[grep("^env.*R", names(asr$vparameters))]
asr$loglik
summary.asreml(asr)$aic
summary.asreml(asr)$bic
sum(diag(Lam %*% t(Lam)))/sum(diag(Lam %*% t(Lam) + Psi)) # variance explained
# True vs estimated
cor(diag(Lam %*% t(Lam) + Psi),diag(De)) # genetic variance
cor(estR, diag(R)) # residual variance 
cor(ug, rowMeans(estInt.reg)) # implicit main effect
mean(diag(cor(matrix(u, ncol = p), estInt.reg + estInt.diag))) # GE effects

# --- Other checks ----
# Check heritability 
mean(diag(cor(matrix(y, ncol = p), matrix(df.MET$u, ncol = p)))^2)
# Compare correlations in original Ce with cors among GE effects 
cor(c(Ce), c(cor(matrix(u, ncol = p))))
cor(c(Ce), c(var(matrix(y, ncol = p))))
range(c(Ce)-c(cor(matrix(u, ncol = p))))
range(c(Ce)-c(cor(matrix(y, ncol = p))))
hist(Ce[upper.tri(Ce)])
hist(cor(matrix(u, ncol = p))[upper.tri(cor(matrix(u, ncol = p)))])
hist(cor(matrix(y, ncol = p))[upper.tri(cor(matrix(y, ncol = p)))])
plot(var(matrix(u, ncol = p)), Ge); abline(a=0, b=1)




########################################################################
# Approach 2: Simulate genotype slopes using AlphaSimR
########################################################################

# --- Use variables simulated in the previous example ---
mu
tau 
X 
Ce
De
Ge
S 
L
e

# --- Simulate slopes with AlphaSimR ---
founderPop = quickHaplo(nInd     = v,
                        nChr     = 1,
                        segSites = 1000,
                        inbred   = TRUE)
SP = SimParam$new(founderPop)

# Simulate k traits as terms with 100 QTLs with mu = 0 and var = 1
SP$addTraitA(nQtlPerChr = 100,
             mean = rep(0, times = k), # set trait mean to 0
             var  = diag(L),           # NOTE: different variance for each term
             corA = diag(k))           # set terms (traits) to be uncorrelated
pop = newPop(founderPop)
slopes = pop@gv
# TO DO: add other non-additive effects

# Add SNP chip
SP$addSnpChip(nSnpPerChr = 900)

# Check traits
apply(slopes, 2, mean)
apply(slopes, 2, var)

# Construct genotype by environment interaction effects
u2 = kronecker(S, diag(v)) %*% c(slopes)

# Obtain genotype main effects
ug2 = rowMeans(matrix(u2, ncol = p))

# Construct an initial MET data frame
df.MET2 = data.frame(env = factor(rep(1:p, each = v * b)),
                     rep = factor(1:b),
                     id  = factor(rep(1:v, each = b)),
                     u   = rep(u2, each = b))
df.MET2 = df.MET2[order(df.MET2$env, df.MET2$rep),]

# --- Assign a randomized complete block design with agricolae ---
# Reorder by trial using design from prevous example
for(i in 1:(p * b)){
  df.MET2[take[i, ], ] = df.MET2[take[i, ], ][design$`1:v`[take[i, ]], ]
}
table(design$`1:v` == df.MET2$id)

# --- Create phenotypes ---
y2 = mu + X %*% tau + df.MET2$u + e
# Check heritability 
mean(diag(cor(matrix(y2, ncol = p), matrix(df.MET2$u, ncol = p)))^2)
# Optional: Assign simulated phenotypes to AlphaSimR pop object
pop@pheno = matrix(y2, ncol = p)

# --- Construct MET dataset ---
df.MET2 = data.frame(env = df.MET2$env,
                     rep = df.MET2$rep,
                     col = df.error$col,
                     row = df.error$row,
                     id  = df.MET2$id,
                     y   = y2,
                     u   = df.MET2$u,
                     e   = e)
head(df.MET2)

# --- Run model ---
# Obtain relationship matrix
Mtt = pullSnpGeno(pop)
Ktt = (Mtt %*% t(Mtt))/10000
diag(Ktt) = diag(Ktt) + 0.00001
Kinv = G.inverse(G = Ktt, sparseform = T)$Ginv
# Model
asr2 = asreml(y ~ 1 + env,
              random    = ~ rr(env, 4):vm(id, Kinv) + diag(env):id,
              residual  = ~ dsum( ~ ar1(col):ar1(row) | env),
              # residual  = ~ dsum( ~ units|env),
              na.action = na.method("include"),
              data = df.MET2)

while (!asr2$converge) {
  asr2 = update.asreml(asr2)
}

# Report
summary(asr2)$varcomp
estInt.reg2 = matrix(asr2$coefficients$random[grep("^rr.*id", 
                rownames(asr2$coefficients$random))][1:(v * p)],ncol = p,byrow = F)
estInt.diag2 = matrix(asr2$coefficients$random[grep("^env.*id", 
                rownames(asr2$coefficients$random))][1:(v * p)],ncol = p,byrow = F)
Lam2 = matrix(asr2$vparameters[grep("rr.*fa", names(asr2$vparameters))], ncol = 4)
Psi2 = diag(asr2$vparameters[grep("^env.*id", names(asr2$vparameters))])
estR2 = asr2$vparameters[grep("^env.*R", names(asr2$vparameters))]
asr2$loglik
summary.asreml(asr2)$aic
summary.asreml(asr2)$bic
sum(diag(Lam2 %*% t(Lam2)))/sum(diag(Lam2 %*% t(Lam2) + Psi2)) # variance explained
# True vs estimated
cor(diag(Lam2 %*% t(Lam2) + Psi2),diag(De)) # genetic variance
cor(estR2, diag(R)) # residual variance 
cor(ug2, rowMeans(estInt.reg2)) # implicit main effect
mean(diag(cor(matrix(u2, ncol = p), estInt.reg2 + estInt.diag2))) # GE effects


# --- Other checks ----
# Check heritability 
mean(diag(cor(matrix(y, ncol = p), matrix(df.MET$u, ncol = p)))^2)
# Compare correlations in original Ce with cors among GE effects 
cor(c(Ce), c(cor(matrix(u, ncol = p))))
cor(c(Ce), c(var(matrix(y, ncol = p))))
range(c(Ce)-c(cor(matrix(u, ncol = p))))
range(c(Ce)-c(cor(matrix(y, ncol = p))))
hist(Ce[upper.tri(Ce)])
hist(cor(matrix(u, ncol = p))[upper.tri(cor(matrix(u, ncol = p)))])
hist(cor(matrix(y, ncol = p))[upper.tri(cor(matrix(y, ncol = p)))])
plot(var(matrix(u, ncol = p)), Ge); abline(a=0, b=1)
