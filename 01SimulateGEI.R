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
# Method 2: structured GxY interaction (TO DO),
# Method 3: structured GxL interaction (TO DO),
# Method 4: structured GxY and GxL interaction (TO DO).
#
# In all methods, GxE interaction is simulated using 
# an extension of Hardin et al. (2013) with the 
# ability to skew the distribution of correlations


# ---- Clean environment and load functions and packages ----

rm(list = ls())
source("00SimulateFunctions.R")

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

# Measures of variance explained
De = diag(1/rgamma(p, shape = 5, rate = 5))
# De = diag(rgamma(n = p, shape = 1.5, scale = 1))
Ge = sqrt(De) %*% Ce %*% sqrt(De)
measureGEI(Ge, prop = T)

# ====================================
#   Measures of variance explained
# ====================================
#                          Proportion
# Main effect variance   : 0.48
# Interaction variance   : 0.52
# ------------------------------------
# Non-crossover variance : 0.51
# Crossover variance     : 0.49
# ====================================


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

# Measures of variance explained
Ge = sqrt(De) %*% Ce %*% sqrt(De)
measureGEI(Ge, prop = T)

# ====================================
#   Measures of variance explained
# ====================================
#                          Proportion
# Main effect variance   : 0.58
# Interaction variance   : 0.42
# ------------------------------------
# Non-crossover variance : 0.64
# Crossover variance     : 0.36
# ====================================


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

# Measures of variance explained
Ge = sqrt(De) %*% Ce %*% sqrt(De)
measureGEI(Ge, prop = T)

# ====================================
#   Measures of variance explained
# ====================================
#                          Proportion
# Main effect variance   : 0.37
# Interaction variance   : 0.63
# ------------------------------------
# Non-crossover variance : 0.42
# Crossover variance     : 0.58
# ====================================


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

# Measures of variance explained
Ge = sqrt(De) %*% Ce %*% sqrt(De)
measureGEI(Ge, prop = T)

# ====================================
#   Measures of variance explained
# ====================================
#                          Proportion
# Main effect variance   : 0.1
# Interaction variance   : 0.9
# ------------------------------------
# Non-crossover variance : 0.11
# Crossover variance     : 0.89
# ====================================


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


# TO DO
# Method 2: structured GxY interaction
# Method 3: structured GxL interaction
# Method 4: structured GxY and GxL interaction