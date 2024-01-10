# Script name: Simulate TPE and generate samples
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# Description:
# This script simulates a TPE dataset with 1000 environments. 
# In addition, lists of 1000 samples of 10, 20 and 50 environments
# are generated with random sampling to create MET datasets for
# analysis in 04ModelComparison script.


# ---- Clean environment and load functions and packages ----

rm(list = ls())
source("../00SimulateFunctions.R")

# Load packages
packages <- c("FieldSimR","asreml","ggpubr","dplyr")
# install.packages(pkgs = packages, dependencies = TRUE)
lapply(packages, library, character.only = TRUE)


# ---- Input parameters and data frame ----

# Scenario parameters
set.seed(123)
scenario    = "ModerateGEI" # Scenario name
subscenario = "noSpatial"   # Sub scenario name
n_genos     = 400           # Number of genotypes
n_reps      = 2             # Number of replicates per environment
n_envs      = 1000          # Number of environments in TPE
n_samples   = 1000          # Number of samples of environments

# Simulate between-environment genetic correlation matrix
source(paste0("simCmat_",scenario,".R"))

#------------------------------------------------------
#-- Simulate the environment genetic variances 
#------------------------------------------------------
# Simulate environment genetic variances
De <- readRDS("De.rds") # use pre-simulated variances
# De <- rgamma(n = n_envs, shape = 1.5, scale = 1)
hist(De, breaks = 100)
# Construct between-environment genetic variance matrix
Ge <- sqrt(diag(De)) %*% Ce %*% sqrt(diag(De))
# Check
mean(Ce[upper.tri(Ce)])
measureGEI(Ge, prop = T) # proportions
measureGEI(Ge, prop = F) # variance
# plotCmat(Ce)$heatmap

#------------------------------------------------------
#-- Simulate phenotypes
#------------------------------------------------------
df <- simPheno(Cmat    = Ce,
               De      = De, 
               n_genos = n_genos, 
               n_reps  = n_reps, 
               mu      = 4,      # overall trait mean
               varE    = 1,      # environment main effect var
               h2      = 0.4,    # plot-level heritability 
               prop_spatial = 0, # if spatial, set prop > 0
               n_cols  = 40,
               n_rows  = 20)

#------------------------------------------------------
#-- Create 1000 samples of 50, 20, 10, and 5 envs
#------------------------------------------------------
# Create 1000 samples of 50 environments from TPE
samples50 <- lapply(X = 1:n_samples, 
                    function(X) sample(1:n_envs, size = 50, replace = F))

# Ensure all subsets of Ge's have positive first loading for 50 environments
resampled <- 0
for(i in 1:length(samples50)){
  cat("Sample:",i,"\n")
  subGe50 <- df$Ge[samples50[[i]],samples50[[i]]]
  subGe50 <- (length(unique(sign(svd(subGe50)$u[,1]))) > 1)
  subGe20 <- df$Ge[samples50[[i]][1:20],samples50[[i]][1:20]]
  subGe20 <- (length(unique(sign(svd(subGe20)$u[,1]))) > 1)
  subGe10 <- df$Ge[samples50[[i]][1:10],samples50[[i]][1:10]]
  subGe10 <- (length(unique(sign(svd(subGe10)$u[,1]))) > 1)
  subGe5 <- df$Ge[samples50[[i]][1:5],samples50[[i]][1:5]]
  subGe5 <- (length(unique(sign(svd(subGe5)$u[,1]))) > 1)
  if (subGe50 == T | subGe20 == T | subGe10 == T | subGe5 == T) {
    cat("Resample:",i,"\n")
    resampled <- resampled + 1
    newSample <- sample(1:n_samples, size = 50, replace = F)
    subGe50 <- df$Ge[newSample,newSample]
    subGe50 <- (length(unique(sign(svd(subGe50)$u[,1]))) > 1)
    subGe20 <- df$Ge[newSample[1:20],newSample[1:20]]
    subGe20 <- (length(unique(sign(svd(subGe20)$u[,1]))) > 1)
    subGe10 <- df$Ge[newSample[1:10],newSample[1:10]]
    subGe10 <- (length(unique(sign(svd(subGe10)$u[,1]))) > 1)
    subGe5 <- df$Ge[newSample[1:5],newSample[1:5]]
    subGe5 <- (length(unique(sign(svd(subGe5)$u[,1]))) > 1)
    while(subGe50 == T | subGe20 == T | subGe10 == T | subGe5 == T) {
      cat("Resample:",i,"\n")
      newSample <- sample(1:n_samples, size = 50, replace = F)
      subGe50 <- df$Ge[newSample,newSample]
      subGe50 <- df$Ge[newSample,newSample]
      subGe50 <- (length(unique(sign(svd(subGe50)$u[,1]))) > 1)
      subGe20 <- df$Ge[newSample[1:20],newSample[1:20]]
      subGe20 <- (length(unique(sign(svd(subGe20)$u[,1]))) > 1)
      subGe10 <- df$Ge[newSample[1:10],newSample[1:10]]
      subGe10 <- (length(unique(sign(svd(subGe10)$u[,1]))) > 1)
      subGe5 <- df$Ge[newSample[1:5],newSample[1:5]]
      subGe5 <- (length(unique(sign(svd(subGe5)$u[,1]))) > 1)
    }
    samples50[[i]] <- newSample
  }
}
resampled
# Subset envs of 20, 10 and 5 from 50 sampled envs
samples20 <- lapply(X = samples50, function(X) X[1:20])
samples10 <- lapply(X = samples50, function(X) X[1:10])
samples5 <- lapply(X = samples50, function(X) X[1:5])
hist(unlist(samples50),breaks = 100)

#------------------------------------------------------
#-- Get summary stats of samples
#------------------------------------------------------
# Output data frame
nSamples = 1000; nScens = 3; nVars = 7
df.smp <- data.frame(matrix(NA, nrow = nSamples*nScens*nVars, ncol = 4))
colnames(df.smp) <- c("Scenario", "Variable", "Sample", "Value")    # GEI variances
df.smp$Variable <- c("varG", "varE", "mainG", "GEI", "nGEI", "cGEI", "h2")
df.smp$Scenario <- c(rep("50 envs",nVars), rep("20 envs",nVars), rep("10 envs",nVars), rep("5 envs",nVars))
df.smp$Sample <- rep(1:nSamples, each = nScens*nVars)
df.smp$Value <- as.numeric(df.smp$Value)

for(i in 1:nSamples){
  cat("Sample:",i,"\n")
  # 50 envs
  subGe <- df$Ge[samples50[[i]],samples50[[i]]]
  df.smp[df.smp$Sample == i & df.smp$Scenario == "50 envs" & df.smp$Variable == "varG",4] <- mean(df$varG[samples50[[i]]])
  df.smp[df.smp$Sample == i & df.smp$Scenario == "50 envs" & df.smp$Variable == "varE",4] <- mean(df$varE[samples50[[i]]])
  df.smp[df.smp$Sample == i & df.smp$Scenario == "50 envs" & df.smp$Variable == "mainG",4] <- measureGEI(subGe,prop = T)$Gvar
  df.smp[df.smp$Sample == i & df.smp$Scenario == "50 envs" & df.smp$Variable == "GEI",4] <- measureGEI(subGe,prop = T)$GEIvar
  df.smp[df.smp$Sample == i & df.smp$Scenario == "50 envs" & df.smp$Variable == "nGEI",4] <- measureGEI(subGe,prop = T,disentangle = T)$nGEI
  df.smp[df.smp$Sample == i & df.smp$Scenario == "50 envs" & df.smp$Variable == "cGEI",4] <- measureGEI(subGe,prop = T,disentangle = T)$cGEI
  df.smp[df.smp$Sample == i & df.smp$Scenario == "50 envs" & df.smp$Variable == "h2",4] <- mean(df$varG[samples50[[i]]]/(df$varG[samples50[[i]]]+df$varE[samples50[[i]]]))
  # 20 envs
  subGe <- df$Ge[samples20[[i]],samples20[[i]]]
  df.smp[df.smp$Sample == i & df.smp$Scenario == "20 envs" & df.smp$Variable == "varG",4] <- mean(df$varG[samples20[[i]]])
  df.smp[df.smp$Sample == i & df.smp$Scenario == "20 envs" & df.smp$Variable == "varE",4] <- mean(df$varE[samples20[[i]]])
  df.smp[df.smp$Sample == i & df.smp$Scenario == "20 envs" & df.smp$Variable == "mainG",4] <- measureGEI(subGe,prop = T)$Gvar
  df.smp[df.smp$Sample == i & df.smp$Scenario == "20 envs" & df.smp$Variable == "GEI",4] <- measureGEI(subGe,prop = T)$GEIvar
  df.smp[df.smp$Sample == i & df.smp$Scenario == "20 envs" & df.smp$Variable == "nGEI",4] <- measureGEI(subGe,prop = T,disentangle = T)$nGEI
  df.smp[df.smp$Sample == i & df.smp$Scenario == "20 envs" & df.smp$Variable == "cGEI",4] <- measureGEI(subGe,prop = T,disentangle = T)$cGEI
  df.smp[df.smp$Sample == i & df.smp$Scenario == "20 envs" & df.smp$Variable == "h2",4] <- mean(df$varG[samples20[[i]]]/(df$varG[samples20[[i]]]+df$varE[samples20[[i]]]))
  # 10 envs
  subGe <- df$Ge[samples10[[i]],samples10[[i]]]
  df.smp[df.smp$Sample == i & df.smp$Scenario == "10 envs" & df.smp$Variable == "varG",4] <- mean(df$varG[samples10[[i]]])
  df.smp[df.smp$Sample == i & df.smp$Scenario == "10 envs" & df.smp$Variable == "varE",4] <- mean(df$varE[samples10[[i]]])
  df.smp[df.smp$Sample == i & df.smp$Scenario == "10 envs" & df.smp$Variable == "mainG",4] <- measureGEI(subGe,prop = T)$Gvar
  df.smp[df.smp$Sample == i & df.smp$Scenario == "10 envs" & df.smp$Variable == "GEI",4] <- measureGEI(subGe,prop = T)$GEIvar
  df.smp[df.smp$Sample == i & df.smp$Scenario == "10 envs" & df.smp$Variable == "nGEI",4] <- measureGEI(subGe,prop = T,disentangle = T)$nGEI
  df.smp[df.smp$Sample == i & df.smp$Scenario == "10 envs" & df.smp$Variable == "cGEI",4] <- measureGEI(subGe,prop = T,disentangle = T)$cGEI
  df.smp[df.smp$Sample == i & df.smp$Scenario == "10 envs" & df.smp$Variable == "h2",4] <- mean(df$varG[samples10[[i]]]/(df$varG[samples10[[i]]]+df$varE[samples10[[i]]]))
  # 5 envs
  subGe <- df$Ge[samples5[[i]],samples5[[i]]]
  df.smp[df.smp$Sample == i & df.smp$Scenario == "5 envs" & df.smp$Variable == "varG",4] <- mean(df$varG[samples5[[i]]])
  df.smp[df.smp$Sample == i & df.smp$Scenario == "5 envs" & df.smp$Variable == "varE",4] <- mean(df$varE[samples5[[i]]])
  df.smp[df.smp$Sample == i & df.smp$Scenario == "5 envs" & df.smp$Variable == "mainG",4] <- measureGEI(subGe,prop = T)$Gvar
  df.smp[df.smp$Sample == i & df.smp$Scenario == "5 envs" & df.smp$Variable == "GEI",4] <- measureGEI(subGe,prop = T)$GEIvar
  df.smp[df.smp$Sample == i & df.smp$Scenario == "5 envs" & df.smp$Variable == "nGEI",4] <- measureGEI(subGe,prop = T,disentangle = T)$nGEI
  df.smp[df.smp$Sample == i & df.smp$Scenario == "5 envs" & df.smp$Variable == "cGEI",4] <- measureGEI(subGe,prop = T,disentangle = T)$cGEI
  df.smp[df.smp$Sample == i & df.smp$Scenario == "5 envs" & df.smp$Variable == "h2",4] <- mean(df$varG[samples5[[i]]]/(df$varG[samples5[[i]]]+df$varE[samples5[[i]]]))
}

#-----------------------------------------------------------------------
# Supplementary plot 1
#-----------------------------------------------------------------------
# Rename and refactorize variables
levels(factor(df.smp$Variable))
df.smp$Variable <- recode_factor(df.smp$Variable, 
                                 `varG` = "varG", 
                                 `varE` = "varE",
                                 `mainG` = "mainG",
                                 `GEI`  = "GEI",
                                 `nGEI` = "nGEI",
                                 `cGEI` = "cGEI",
                                 `h2`   = "h2")
# Rename and refactorize variables
levels(factor(df.smp$Scenario))
df.smp$Scenario <- recode_factor(df.smp$Scenario, `5 envs` = "5 envs", `10 envs` = "10 envs",
                                 `20 envs` = "20 envs",`50 envs` = "50 envs")

# Add addtional points for true values
df.true <- data.frame("Scenario" = c(rep("50 envs",nVars), rep("20 envs",nVars), rep("10 envs",nVars), rep("5 envs",nVars)), 
                      "Variable" = c("varG", "varE", "mainG", "GEI", "nGEI", "cGEI", "h2"),
                      "Value" = as.numeric(NA))
df.true$Value <- c(mean(df$varG),mean(df$varE),
                   measureGEI(df$Ge,prop = T)$Gvar,
                   measureGEI(df$Ge,prop = T)$GEIvar,
                   # NA,NA,
                   measureGEI(df$Ge,prop = T,disentangle = T)$nGEI, #<<<<
                   measureGEI(df$Ge,prop = T,disentangle = T)$cGEI, #<<<<
                   mean(df$varG/(df$varG+df$varE)))
# Rename and refactorize variables
df.true$Scenario <- recode_factor(df.true$Scenario, `5 envs` = "5 envs", `10 envs` = "10 envs",
                                 `20 envs` = "20 envs",`50 envs` = "50 envs")

# Plot 1
temp <- df.true %>%
  filter((Variable %in% c("mainG", "GEI", "nGEI", "cGEI", "h2"))) %>%
  droplevels()
(a <- df.smp %>%
    filter((Variable %in% c("mainG", "GEI", "nGEI", "cGEI", "h2"))) %>%
    droplevels() %>% 
    ggplot(aes(x = Variable, y = Value, fill = Variable)) +
    facet_grid(~Scenario, scales = "free") + 
    geom_boxplot() +
    # geom_point(temp,aes(x = Variable, y = Value, fill = Variable),shape=20, size=3, color="red") + 
    geom_point(data = temp, aes(x = Variable, y = Value), shape = 4, size = 3) +
    # stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="black", fill="black") +
    ylim(0,1) +
    theme_bw(base_size = 16, base_family = "sans") +
    theme(panel.background = element_blank(),
          # legend = element_blank(),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 14),
          legend.position = "none",
          plot.title = element_text(size = 15, face = "bold"),
          axis.text  = element_text(size = 14,colour = "black"),
          axis.title = element_text(size = 18),
          # axis.title.x = element_blank(),
          axis.text.x  = element_text(vjust = 0.7, hjust = 0.5, angle = 45),
          strip.text   = element_text(face = "bold", size = 18, colour = 'black'))
)

# Plot 2
temp <- df.true %>%
  filter((Variable %in% c("varG", "varE"))) %>%
  droplevels()
(b <- df.smp %>%
    filter((Variable %in% c("varG", "varE"))) %>%
    droplevels() %>% 
    ggplot(aes(x = Variable, y = Value, fill = Variable)) +
    facet_grid(~Scenario, scales = "free") + 
    geom_boxplot() +
    # geom_point(temp,aes(x = Variable, y = Value, fill = Variable),shape=20, size=3, color="red") + 
    geom_point(data = temp, aes(x = Variable, y = Value), shape = 4, size = 3) +
    # stat_summary(fun.y=mean, geom="point", shape=20, size=3, color="black", fill="black") +
    ylab("Variance") +
    ylim(0,20) +
    theme_bw(base_size = 16, base_family = "sans") +
    theme(panel.background = element_blank(),
          # legend = element_blank(),
          legend.title = element_text(size = 14, face = "bold"),
          legend.text  = element_text(size = 14),
          legend.position = "none",
          plot.title = element_text(size = 15, face = "bold"),
          axis.text  = element_text(size = 14,colour = "black"),
          axis.title = element_text(size = 18),
          # axis.title.x = element_blank(),
          axis.text.x  = element_text(vjust = 0.7, hjust = 0.5, angle = 45),
          strip.text   = element_text(face = "bold", size = 18, colour = 'black'))
)

(p <- ggarrange(a, b, nrow = 2, legend = "none", align = "hv"))
ggsave(plot = p, filename = paste0("SupPlot-samples-",scenario,".pdf"), 
       width = 6, height = 6, scale = 1.7)

####################################################################################
# Save data frame object - contains true values, full TPE and lists of 1000 samples
# of 10, 20 and 50 environments
df[["samples5"]]  <- samples5
df[["samples10"]] <- samples10
df[["samples20"]] <- samples20
df[["samples50"]] <- samples50
save(list = "df", file = paste0("samples",scenario,".RData"))

