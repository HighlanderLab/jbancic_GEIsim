# Script name: Basic framework to simulate phenotypic data (Appendix A)
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# Description:
# This script demonstrates simulation of phenotypic data to construct
# a MET dataset. This code is also provided in Appendix 1.

# ---- Clean environment and load functions and packages ----
rm(list = ls())
source("00SimulateFunctions.R")

# Load packages
packages <- c("AlphaSimR", "FieldSimR", "agricolae", 
              "asreml", "ASRgenomics", "ggplot2")
# install.packages(pkgs = packages, dependencies = TRUE)
lapply(packages, library, character.only = TRUE)


#> Initial parameters
set.seed(123) # do not change
p  = 10   # Environments
b  = 2    # Blocks
v  = 200  # Genotypes
mu = 4    # Trait mean
H2 = 0.3  # Plot-level heritability
k  = 7    # No. of multiplicative terms

#> Simulate environmental effects
X   = kronecker(diag(p), rep(1, v * b))
tau = rnorm(p, 0, 1)

#> Simulate GE effects
# 1. Simulate Ge
Ce = matrix(0, p, p)
Ce[upper.tri(Ce,F)] = runif(p * (p - 1)/2, 0.4, 1)
Ce = Ce + t(Ce); diag(Ce) = 1
De = diag(1/rgamma(p, shape = 5, rate = 5))
Ge = sqrt(De) %*% Ce %*% sqrt(De)

# 2. Decompose Ge and take first k terms
U = svd(Ge)$u[,1:k]
L = diag(svd(Ge)$d[1:k])

# 3. Obtain covariates and slopes
S = U
slopes = scale(matrix(rnorm(k * v), ncol = k))
slopes = c(matrix(slopes, ncol = k) %*% sqrt(L))

# 4. Construct GE effects 
Z = kronecker(diag(v * p), rep(1, b))
u = kronecker(S, diag(v)) %*% slopes

# Obtain genotype main effects
ug = rowMeans(matrix(u, ncol = p))

#> Simulate plot errors
H2 = abs(rnorm(p, H2, 0.1))
H2[H2 < 0] = 0; H2[H2 > 1] = 1
R  = diag(diag(De) / H2 - diag(De))
e  = c(matrix(rnorm(p*b*v), ncol = p) %*% sqrt(R))

#> Create phenotypes
y = mu + X %*% tau + Z %*% u + e

# Construct MET dataset
df.MET = data.frame(
  env = factor(rep(1:p, each = v * b)),
  rep = factor(1:b),
  id  = factor(rep(1:v, each = b)),
  y = y,
  u = rep(u, each = b),
  e = e
  )

# Order MET and randomise by trial
df.MET = df.MET[order(df.MET$env, df.MET$rep), ]
for(i in 0:(p * b - 1)){
  df.MET[i * v + 1:v,] = df.MET[i * v + 1:v, ][sample(1:v, v, F), ]
}

#> Run model
asr = asreml(y ~ 1 + env,
             random    = ~ rr(env, 4):id + 
               diag(env):id,
             residual  = ~ dsum(~units|env),
             na.action = na.method("include"),
             data      = df.MET)

while (!asr$converge) {
  asr = update.asreml(asr)
}
# save.image("Paper_ToyExample.RData")

########################################################################
# Report model results
summary(asr)$varcomp
estInt.reg  <- matrix(asr$coefficients$random[grep("^rr.*id", rownames(asr$coefficients$random))][1:(v*p)],ncol = p,byrow = F)
estInt.diag <- matrix(asr$coefficients$random[grep("^env.*id", rownames(asr$coefficients$random))][1:(v*p)],ncol = p,byrow = F)
Lam <- matrix(asr$vparameters[grep("rr.*fa", names(asr$vparameters))], ncol = 4)
Psi <- diag(asr$vparameters[grep("^env.*id", names(asr$vparameters))])
estR <- asr$vparameters[grep("^env.*R", names(asr$vparameters))]
asr$loglik
summary.asreml(asr)$aic
summary.asreml(asr)$bic
sum(diag(Lam %*% t(Lam)))/sum(diag(Lam %*% t(Lam) + Psi)) # VAF
cor(diag(Lam %*% t(Lam) + Psi),diag(De)) # genetic var
cor(estR,diag(R)) # residual var
sqrt(sum((Ge[upper.tri(Ge,diag = T)] - (Lam %*% t(Lam) + Psi)[upper.tri(Ge,diag = T)])^2)/length(Ge[upper.tri(Ge,F)])) # RMSD Ge cors
sqrt(sum((Ce[upper.tri(Ge,diag = T)] - cov2cor(Lam %*% t(Lam) + Psi)[upper.tri(Ge,diag = T)])^2)/length(Ge[upper.tri(Ge,F)])) # RMSD Ce cors
cor(ug,rowMeans(estInt.reg)) # implicit main effect
cor(ug,rowMeans(estInt.reg + estInt.diag))
mean(diag(cor(matrix(u,ncol=p),estInt.reg)))
mean(diag(cor(matrix(u,ncol=p),estInt.reg + estInt.diag)))

########################################################################
## Check different components
# Check variance matrix
range(diag(k) - t(S) %*% S)
range(cov2cor(Ge) - cov2cor(S %*% L %*% t(S)))
# Check genotype slopes
apply(matrix(slopes, ncol = p), 2, mean)
apply(matrix(slopes, ncol = p), 2, var)
plot(diag(var(matrix(slopes, ncol = k))), diag(L)); abline(a=0, b=1)
# Check MET
tapply(df.MET$y,df.MET$env,mean); mean(df.MET$y)
tapply(df.MET$u,df.MET$env,var);  mean(tapply(df.MET$u,df.MET$env,var))
tapply(df.MET$e,df.MET$env,var);  mean(tapply(df.MET$e,df.MET$env,var))
# Check plot-level heritability
tapply(df.MET$u,df.MET$env,var)/(tapply(df.MET$u,df.MET$env,var)+tapply(df.MET$e,df.MET$env,var))
mean(tapply(df.MET$u,df.MET$env,var)/(tapply(df.MET$u,df.MET$env,var)+tapply(df.MET$e,df.MET$env,var)))

########################################################################
## Plot genetic matrices
plotCmat(Ce)$heat
plotCmat(Ge[plotCmat(Ce)$order,plotCmat(Ce)$order])$heat

########################################################################
library(ggplot2)
## OP-RMSD plot 
OP   <- mean(S[,1]) * slopes[1:v] # Lam * scores
dev  <- matrix(slopes, ncol = k)[,2:k] %*% t(S[,2:k]) # higher order factors
RMSD <- sqrt(rowMeans(dev^2)) # square root of means
df   <- data.frame(geno = as.factor(1:v), OP, RMSD)

(p.OP_RMSD <- ggplot(df, aes(x=RMSD, y=OP)) +
    geom_point(size = 1.5, colour = "gray40") +
    geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.7) +
    labs(y = "Overall performance",x = "Stability"))

########################################################################
## Regression plot 
temp <- data.frame(Env    = rep(1:p, each = v),
                   geno   = as.factor(1:v),
                   Lam    = rep(S[,1:k], each = v),
                   f      = c(rep(1,p) %x% matrix(slopes, ncol = k)[,1:k]),
                   Source = as.factor(rep(paste("Term",1:k, sep = " "), each = v*p)),
                   CVE    = c(matrix(slopes, ncol = k)[,2:k] %*% t(S[,2:k])))
# Rename factor names
temp$Source <- factor(temp$Source, levels = paste("Term",1:p, sep = " "))
levels(temp$Source)[7] <- "Term k"
take <- c(paste("Term",c(1:2), sep = " "),"Term k")
# Set up horizontal and vertical lines
meanLoad <- -mean(S[,1])
vlines <- data.frame(Source = unique(temp$Source)[c(1:2,7)], 
                     vline  = c(meanLoad,0,0))
hlines <- data.frame(Source = unique(temp$Source)[c(1:2,7)], 
                     hline  = c(0,0,0))

(p.reg <- ggplot(temp[temp$Source %in% take,], aes(y = Lam*f, x = Lam)) +
    geom_line(aes(y= Lam*f, x = -Lam, group = geno)) +
    geom_vline(data = vlines, 
               aes(xintercept = vline), color = "black", linetype = "dotted", size = 0.7) +
    geom_hline(data = hlines, 
               aes(yintercept = hline), color = "black", linetype = "dotted", size = 0.7) +
    facet_wrap(~ Source, scales = "free_x",ncol = 4) +
    labs(x = "Environmental covariates",
         y = "True genetic effects"))
# ggsave(plot = p.reg, filename ="regPlot.jpeg", width = 7, height = 3, scale = 2)
