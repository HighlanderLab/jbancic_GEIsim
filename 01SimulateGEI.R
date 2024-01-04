# Script name: A framework to simulate GxE interaction
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# Description:
# This script demonstrates simulation of GxE interaction.
# Four methods are used to simulate GxE interaction:
#
# Method 1: structured GxE interaction,
# Method 2: structured GxY interaction,
# Method 3: structured GxL interaction,
# Method 4: structured GxY and GxL interaction.
#
# In all methods, GxE interaction is simulated using 
# an extension of Hardin et al. (2013) with the 
# ability to skew the distribution of correlations


# ---- Clean environment and load functions and packages ----

rm(list = ls())
source("00SimulateFunctions.R")
# install.packages(pkgs = "AlphaSimR")
# library(package = "AlphaSimR")

# ---- Specify parameters ----
p  = 100  # Number of environments


#######################################################
# Method 1: structured GxE interaction
#######################################################

# We consider 3 examples
# 1. A simple example
# 2. Three examples with skewed correlations
#      i) Low crossover GxE interaction
#     ii) Moderate crossover GxE interaction
#    iii) High crossover GxE interaction
# 3. An example with multiple environment groups


# ---- 1. Simple example ----
# stop function: out of range, respecify epsilon and/or mean_cor
# is abs(supplied) value > 0.99 - abs(mean_cor)
Ce = simCmat(n_envs = p,
               mean_cor = 0.999-0.5,
               rank = 6, 
               epsilon = 0.5)
plotCmat(Ce, den_order = T)$heat
plotCmat(Ce, den_order = T)$hist
plotCmat(Ce, den_order = T)$data
plotCmat(Ce, den_order = T)$order

# Calculate measures of variance explained
De = diag(rgamma(n = p, shape = 1.5, scale = 1))
Ge = sqrt(De) %*% Ce %*% sqrt(De)
plot(diag(Ge), svd(Ge)$u[,1]^2)
plot(svd(Ce)$u[,1], svd(Ge)$u[,1])
# G
mean(Ce)
100*sum(svd(Ce)$d*colMeans(svd(Ce)$u)^2)/sum(diag(Ce))
measureGEI(Ce, prop = F)$Gvar
# ~70% w/ homo gen vars
measureGEI(Ge)$Gvar
# ~60% G w/ het gen vars
# nGEI
measureGEI(Ce, disentangle = T)$nGEI
# ~70% nGEI w/ homo gen vars, so nothing that is marginal
measureGEI(Ge, disentangle = T)$nGEI
# ~70% nGEI w/ het gen vars, so ~10% that is marginal
#cGEI
measureGEI(Ce, disentangle = T)$cGEI
# ~30% cGEI w/ homo gen vars
measureGEI(Ge, disentangle = T)$cGEI
# ~30% cGEI w/ het gen vars
# P = diag(100) - matrix(1, ncol = 100, nrow = 100)/100 # av-semi variance
# sum(diag(Ce %*% P)) / 99 # 0.300772
# sum(diag(matrix(mean(Ce), ncol = 100, nrow = 100) %*% P)) / 99 # 0



# ---- 2. Example with skewed distribution of correlations ----
# i) Low crossover GEI
Ce = simCmat(n_envs = p,
               mean_cor = 0.5,
               epsilon = 0.999-0.5,
               rank = 6,
               skew = 0.5)
plotCmat(Ce, den_order = T)$heatmap
plotCmat(Ce, den_order = T)$hist
plotCmat(Ce, den_order = T)$data
plotCmat(Ce, den_order = T)$order

# measures of variance explained
Ge = sqrt(De) %*% Ce %*% sqrt(De)
# G
100*sum(svd(Ce)$d*colMeans(svd(Ce)$u)^2)/sum(diag(Ce))
measureGEI(Ce)$Gvar
# ~75% G w/ homo gen vars
measureGEI(Ge)$Gvar
# ~65% G w/ het gen vars
# nGEI
measureGEI(Ce, disentangle = T)$nGEI
# ~75% nGEI w/ homo gen vars, so nothing that is marginal
measureGEI(Ge, disentangle = T)$nGEI
# ~75% nGEI w/ het gen vars, so ~10% that is marginal
#cGEI
measureGEI(Ce, disentangle = T)$cGEI
# ~25% cGEI w/ homo gen vars
measureGEI(Ge, disentangle = T)$cGEI
# ~25% cGEI w/ het gen vars


# ii) Moderate crossover GxE
Ce = simCmat(n_envs = p,
               mean_cor = 0.2,
               epsilon = 0.95-0.2,
               rank = 6,
               skew = 0.5)
plotCmat(Ce, den_order = T)$heatmap
plotCmat(Ce, den_order = T)$hist
plotCmat(Ce, den_order = T)$data
plotCmat(Ce, den_order = T)$order

# measures of variance explained
Ge = sqrt(De) %*% Ce %*% sqrt(De)
# G
100*sum(svd(Ce)$d*colMeans(svd(Ce)$u)^2)/sum(diag(Ce))
measureGEI(Ce)$Gvar
# ~60% G w/ homo gen vars
measureGEI(Ge)$Gvar
# ~50% G w/ het gen vars
# nGEI
measureGEI(Ce, disentangle = T)$nGEI
# ~60% nGEI w/ homo gen vars, so practically nothing that is marginal
measureGEI(Ge, disentangle = T)$nGEI
# ~60% nGEI w/ het gen vars, so ~10% that is marginal
#cGEI
measureGEI(Ce, disentangle = T)$cGEI
# ~40% cGEI w/ homo gen vars
measureGEI(Ge, disentangle = T)$cGEI
# ~40% cGEI w/ het gen vars


# iii) High crossover GxE
Ce = simCmat(n_envs = p,
               mean_cor = 0,
               epsilon = 1,
               rank = 6,
               skew = 0.35)
plotCmat(Ce, den_order = T)$heatmap
plotCmat(Ce, den_order = T)$hist
plotCmat(Ce, den_order = T)$data
plotCmat(Ce, den_order = T)$order

# measures of variance explained
Ge = sqrt(De) %*% Ce %*% sqrt(De)
# G
100*sum(svd(Ce)$d*colMeans(svd(Ce)$u)^2)/sum(diag(Ce))
measureGEI(Ce)$Gvar
# ~40% G w/ homo gen vars
measureGEI(Ge)$Gvar
# ~35% G w/ het gen vars
# nGEI
measureGEI(Ce, disentangle = T)$nGEI
# ~40% nGEI w/ homo gen vars, so practically nothing that is marginal
measureGEI(Ge, disentangle = T)$nGEI
# ~40% nGEI w/ het gen vars, so ~5% that is marginal
#cGEI
measureGEI(Ce, disentangle = T)$cGEI
# ~60% cGEI w/ homo gen vars
measureGEI(Ge, disentangle = T)$cGEI
# ~60% cGEI w/ het gen vars


# ---- 3. Example with multiple environment groups ----
Ce = simCmat(n_groups = 2,
                n_envs = c(50,150),
                mean_cor = 0.5,
                delta = 0.3,
                epsilon = 0.4,
                skew = 0.5)
plotCmat(cor_mat = Ce, den_order = T, groups = list(1:50,51:200))$heatmap
plotCmat(Ce, den_order = T, groups = list(1:50,51:200))$hist
plotCmat(Ce, den_order = T, groups = list(1:50,51:200))$data
plotCmat(Ce, den_order = T, groups = list(1:50,51:200))$order
# functionality to have a different range within and between groups


#######################################################
# Method 2: structured GxY interaction
#######################################################



#######################################################
# Method 3: structured GxL interaction
#######################################################



#######################################################
# Method 4: structured GxY and GxL interaction
#######################################################
