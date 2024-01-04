#!/exports/applications/apps/community/roslin/R/4.2.1-asreml/bin/Rscript
#$ -N noSpH-10
#$ -cwd
#$ -R y
#$ -pe sharedmem 1
#$ -l h_vmem=30G
#$ -l h_rt=50:00:00
#$ -j y
#$ -V
#$ -P roslin_HighlanderLab
#$ -M jbancic@exseed.ed.ac.uk
#$ -m eas

# Script name: Model comparison
#
# Authors: Jon Bancic, Gregor Gorjanc, Daniel Tolhurst
#
# Date Created: 2024-01-01
#
# Description:
# This script demonstrates a single replicate of model 
# comparison example shown in the manuscript.


# ---- Clean environment and load functions and packages ----

rm(list = ls())
source("../00SimulateFunctions.R")

# Load packages
packages <- c("FieldSimR","asreml")
# install.packages(pkgs = packages, dependencies = TRUE)
lapply(packages, library, character.only = TRUE)


# ---- Input parameters and data frame ----

# Scenario parameters
Rep <- 1                   # Replicate number
scenario <- "ModerateGEI"  # Scenario name
subscenario <- "noSpatial" # Sub scenario name
n_genos <- 400             # Number of genotypes in MET dataset
n_reps  <- 2               # Number of replicates per environment
n_envs  <- 10              # Number of environments in MET dataset
n_totenvs <- 1000          # Number of environments in TPE

# Output data frame
df.output <- data.frame(matrix(NA, nrow = 11, ncol = 22))
colnames(df.output) <- c("Scenario", "subScenario", "Envs", "Rep", "Model", 
                         "LogLik", "AIC", "BIC", "runTime", "VAF",   # model criteria
                         "Acc-varG", "Acc-totGe", "Acc-totCe",       # estimates correlations
                         "Acc-0","Acc-1", "Acc-2", "Acc-3", "Acc-4", # accuracies
                         "varG", "varGEI", "nGEI", "cGEI")           # GEI variances
df.output$Model <- c("Maximum", "Expected", "Pheno", "Main", 
                     "Comp", "MDiag", "Diag", "FA1", "FA2", "FA3", "FA4")
df.output$Rep  <- Rep
df.output$Envs <- n_envs
df.output$Scenario  <- scenario
df.output$subScenario <- subscenario

# Run simulation
cat("Working on ", scenario, " rep ", Rep,"\n",sep = "")

#------------------------------------------------------
#-- Load pre-simulated data for MET sample
#------------------------------------------------------
load(paste0("samples",scenario,".RData"))
sample <- df[[paste0("samples",n_envs)]][[Rep]][order(df[[paste0("samples",n_envs)]][[Rep]])]
df.sub <- list(Ce     = df$Ce[sample,sample],
               Ge     = df$Ge[sample,sample],
               trueG  = rowMeans(df$trueGE[,sample]),
               trueGE = df$trueGE[,sample],
               varG   = df$varG[sample],
               varE   = df$varE[sample],
               met.df = droplevels(df$met.df[df$met.df$env %in% sample,]))

#------------------------------------------------------
#-- Simulate phenotypes
#------------------------------------------------------
# Report variance explained by different components in proportion
df.output[df.output$Rep == Rep,]$varG   <- measureGEI(df.sub$Ge, prop = T)$Gvar   # varG
df.output[df.output$Rep == Rep,]$varGEI <- measureGEI(df.sub$Ge, prop = T)$GEIvar # varGEI
df.output[df.output$Rep == Rep,]$nGEI   <- measureGEI(df.sub$Ge, disentangle = T, prop = T)$nGEI # nGEI
df.output[df.output$Rep == Rep,]$cGEI   <- measureGEI(df.sub$Ge, disentangle = T, prop = T)$cGEI # cGEI

# Report true variance components
varG   <- measureGEI(df$Ge,prop = F)$Gvar
varGEI <- measureGEI(df$Ge,prop = F)$GEIvar
varE   <- mean(df$varE)

# Report maximum and expected accuracies
df.output[df.output$Rep == Rep & df.output$Model == "Maximum",]$`Acc-1`  <- cor(df.sub$trueG, df$trueG) # max main effect; equal to the cor between sub and TPE
df.output[df.output$Rep == Rep & df.output$Model == "Expected",]$`Acc-1` <- sqrt(varG/(varG + varGEI/n_envs + varE/n_envs/n_reps)) # expected main eff -- ratio of means
df.output[df.output$Rep == Rep & df.output$Model == "Expected",]$`Acc-2` <- sqrt(varG/(varG + varGEI/n_envs + varE/n_envs/n_reps)) # expected main eff -- ratio of means
df.output[df.output$Rep == Rep & df.output$Model == "Expected",]$`VAF`   <- varG/(varG+varGEI)

#------------------------------------------------------
#-- Run models
#------------------------------------------------------

# NULL: Phenotypic ----------------------
estG   <- with(df.sub$met.df, tapply(y, id, mean))
estGE  <- with(df.sub$met.df, tapply(y, list(id, env), mean))
# estInt <- estGE - matrix(estG, ncol = n_envs, nrow = n_genos)

df.output[df.output$Rep == Rep & df.output$Model == "Pheno",]$`Acc-0` <- cor(df.sub$trueG,estG) 
df.output[df.output$Rep == Rep & df.output$Model == "Pheno",]$`Acc-1` <- cor(df$trueG,estG)
df.output[df.output$Rep == Rep & df.output$Model == "Pheno",]$`Acc-2` <- cor(df$trueG,estG)
df.output[df.output$Rep == Rep & df.output$Model == "Pheno",]$`Acc-3` <- NA
df.output[df.output$Rep == Rep & df.output$Model == "Pheno",]$`Acc-4` <- mean(diag(cor(df.sub$trueGE,estGE)))

# 0. Diagonal model ----------------------
cat("--> Diagonal model \n")
start_time <- Sys.time()
asr.diag <- asreml(y ~ 1 + env,
                   random = ~ diag(env):id,
                   residual = ~dsum(~units|env),
                   na.action = na.method("include"),
                   data = df.sub$met.df)
while (asr.diag$converge != TRUE) {
  asr.diag = update.asreml(asr.diag)
}
end_time <- Sys.time()

# Get effects
# interaction
estInt = asr.diag$coefficients$random[
  grep(pattern = ".*env.*id", rownames(asr.diag$coefficients$random)),]
estInt = matrix(estInt,ncol = n_envs,byrow = F)
# variances
estvarG <- asr.diag$vparameters[grep("id.*env",names(asr.diag$vparameters))]
estGe <- diag(estvarG)

# Report
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$LogLik    <- asr.diag$loglik
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$AIC       <- summary.asreml(asr.diag)$aic
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$BIC       <- summary.asreml(asr.diag)$bic
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$`runTime` <- difftime(end_time, start_time, units = "secs")
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$`Acc-varG`  <- sqrt(sum((diag(df.sub$Ge) - diag(estGe))^2)/length(diag(estGe)))
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$`Acc-totGe` <- sqrt(sum((df.sub$Ge[upper.tri(df.sub$Ge,F)] - estGe[upper.tri(estGe,F)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$`Acc-totCe` <- sqrt(sum((df.sub$Ce[upper.tri(df.sub$Ce,F)] - cov2cor(estGe)[upper.tri(estGe,F)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$`Acc-0` <- NA # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$`Acc-1` <- NA # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$`Acc-2` <- NA
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$`Acc-3` <- NA
df.output[df.output$Rep == Rep & df.output$Model == "Diag",]$`Acc-4` <- mean(diag(cor(matrix(df.sub$trueGE, ncol = n_envs),estInt)))

# 1. Main effects ----------------------
# Model
cat("--> Main effects model \n")
start_time <- Sys.time()
asr.main <- asreml(y ~ 1 + env,
                   random = ~ id,
                   residual = ~dsum(~units|env),
                   na.action = na.method("include"),
                   data = df.sub$met.df)
while (asr.main$converge != TRUE) {
  asr.main = update.asreml(asr.main)
}
end_time <- Sys.time()

# Get effects
estG = asr.main$coefficients$random[
  grep(pattern = ".*id", rownames(asr.main$coefficients$random)),]
names(estG) = sub(pattern = ".*_","", names(estG))
# Report
df.output[df.output$Rep == Rep & df.output$Model == "Main",]$LogLik  <- asr.main$loglik
df.output[df.output$Rep == Rep & df.output$Model == "Main",]$AIC     <- summary.asreml(asr.main)$aic
df.output[df.output$Rep == Rep & df.output$Model == "Main",]$BIC     <- summary.asreml(asr.main)$bic
df.output[df.output$Rep == Rep & df.output$Model == "Main",]$`runTime` <- difftime(end_time, start_time, units = "secs")
df.output[df.output$Rep == Rep & df.output$Model == "Main",]$`Acc-0` <- cor(df.sub$trueG,estG) # explicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "Main",]$`Acc-1` <- cor(df$trueG,estG) # explicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "Main",]$`Acc-2` <- NA
df.output[df.output$Rep == Rep & df.output$Model == "Main",]$`Acc-3` <- NA
df.output[df.output$Rep == Rep & df.output$Model == "Main",]$`Acc-4` <- NA

# 2. Compound Symmetry -----------------
# Model
cat("--> Compound symmetry model \n")
start_time <- Sys.time()
asr.comp <- asreml(y ~ 1 + env,
                   random = ~ id + 
                     idv(env):id,
                   residual = ~dsum(~units|env),
                   na.action = na.method("include"),
                   data = df.sub$met.df)
while (asr.comp$converge != TRUE) {
  asr.comp = update.asreml(asr.comp)
}
end_time <- Sys.time()

# Get effects
# main effect
temp = asr.comp$coefficients$random[
  grep(pattern = ".*id", rownames(asr.comp$coefficients$random)),]
estG = temp[grep(pattern = "env", names(temp),invert = T)]
names(estG) = sub(pattern = ".*_","", names(estG))
# interaction
estInt = asr.comp$coefficients$random[
  grep(pattern = ".*env.*id", rownames(asr.comp$coefficients$random)),]
estInt = matrix(estInt,ncol = n_envs,byrow = F)
# variances
estMainVar <- asr.comp$vparameters[grep("id",names(asr.comp$vparameters))][1]
estvarG <- asr.comp$vparameters[grep("id.*env",names(asr.comp$vparameters))]
estGe <- matrix(estMainVar,ncol = n_envs,nrow = n_envs) + diag(rep(estvarG,times = n_envs))

# Report
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$LogLik  <- asr.comp$loglik
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$AIC     <- summary.asreml(asr.comp)$aic
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$BIC     <- summary.asreml(asr.comp)$bic
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$`runTime` <- difftime(end_time, start_time, units = "secs")
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$VAF     <- estMainVar/(estMainVar+estvarG)
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$`Acc-varG`  <- sqrt(sum((diag(df.sub$Ge) - diag(estGe))^2)/length(diag(estGe)))
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$`Acc-totGe` <- sqrt(sum((df.sub$Ge[upper.tri(df.sub$Ge,F)] - estGe[upper.tri(estGe,F)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$`Acc-totCe` <- sqrt(sum((df.sub$Ce[upper.tri(df.sub$Ce,F)] - cov2cor(estGe)[upper.tri(estGe,F)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$`Acc-0` <- cor(df.sub$trueG,estG) # explicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$`Acc-1` <- cor(df$trueG,estG) # explicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$`Acc-2` <- cor(df$trueG,estG+rowMeans(estInt))
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$`Acc-3` <- mean(cor(df.sub$trueGE,estG))
df.output[df.output$Rep == Rep & df.output$Model == "Comp",]$`Acc-4` <- mean(diag(cor(df.sub$trueGE,(estG+estInt))))

# 3. Main + diagonal model -------------
# Model
cat("--> Main + diag model \n")
start_time <- Sys.time()
asr.md <- asreml(y ~ 1 + env,
                 random = ~ id + 
                   diag(env):id,
                 residual = ~dsum(~units|env),
                 na.action = na.method("include"),
                 data = df.sub$met.df)
while (asr.md$converge != TRUE) {
  asr.md = update.asreml(asr.md)
}
end_time <- Sys.time()

# Get effects
# main effect
temp = asr.md$coefficients$random[
  grep(pattern = ".*id", rownames(asr.md$coefficients$random)),]
estG = temp[grep(pattern = "env", names(temp),invert = T)]
names(estG) = sub(pattern = ".*_","", names(estG))
# interaction
estInt = asr.md$coefficients$random[
  grep(pattern = ".*env.*id", rownames(asr.md$coefficients$random)),]
estInt = matrix(estInt,ncol = n_envs,byrow = F)
# variances
estMainVar <- asr.md$vparameters[grep("id",names(asr.md$vparameters))][1]
estvarG <- asr.md$vparameters[grep("id.*env",names(asr.md$vparameters))]
estGe <- matrix(estMainVar,ncol = n_envs,nrow = n_envs)+diag(estvarG)

# Report
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$LogLik  <- asr.md$loglik
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$AIC     <- summary.asreml(asr.md)$aic
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$BIC     <- summary.asreml(asr.md)$bic
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$`runTime` <- difftime(end_time, start_time, units = "secs")
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$VAF     <- mean(estMainVar/(estMainVar+estvarG))
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$`Acc-varG`  <- sqrt(sum((diag(df.sub$Ge) - diag(estGe))^2)/length(diag(estGe)))
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$`Acc-totGe` <- sqrt(sum((df.sub$Ge[upper.tri(df.sub$Ge,F)] - estGe[upper.tri(estGe,F)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$`Acc-totCe` <- sqrt(sum((df.sub$Ce[upper.tri(df.sub$Ce,F)] - cov2cor(estGe)[upper.tri(estGe,F)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$`Acc-0` <- cor(df.sub$trueG,estG) # explicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$`Acc-1` <- cor(df$trueG,estG) # explicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$`Acc-2` <- cor(df$trueG,estG+rowMeans(estInt))
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$`Acc-3` <- mean(cor(df.sub$trueGE,estG))
df.output[df.output$Rep == Rep & df.output$Model == "MDiag",]$`Acc-4` <- mean(diag(cor(df.sub$trueGE,(estG+estInt))))

# 4. Factor analytic 1 -----------------
# Model
cat("--> FA1 model \n")
start_time <- Sys.time()
asr.fa1 <- asreml(y ~ 1 + env,
                  random = ~ rr(env, 1):id + 
                    diag(env):id,
                  residual = ~dsum(~units|env),
                  na.action = na.method("include"),
                  data = df.sub$met.df,
                  workspace = "1Gb")
while (asr.fa1$converge != TRUE) {
  asr.fa1 = update.asreml(asr.fa1)
}
end_time <- Sys.time()

# Get effects
estInt.reg  <- matrix(asr.fa1$coefficients$random[grep("^rr.*id", rownames(asr.fa1$coefficients$random))][1:(n_genos*n_envs)],ncol = n_envs,byrow = F)
estInt.diag <- matrix(asr.fa1$coefficients$random[grep("^env.*id", rownames(asr.fa1$coefficients$random))][1:(n_genos*n_envs)],ncol = n_envs,byrow = F)
Lam <- matrix(asr.fa1$vparameters[grep("rr.*fa", names(asr.fa1$vparameters))], ncol = 1)
Psi <- diag(asr.fa1$vparameters[grep("^env.*id", names(asr.fa1$vparameters))])
# Report
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$LogLik  <- asr.fa1$loglik
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$AIC     <- summary.asreml(asr.fa1)$aic
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$BIC     <- summary.asreml(asr.fa1)$bic
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$`runTime` <- difftime(end_time, start_time, units = "secs")
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$VAF     <- sum(diag(Lam %*% t(Lam)))/sum(diag(Lam %*% t(Lam) + Psi))
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$`Acc-varG`  <- sqrt(sum((diag(df.sub$Ge) - diag(Lam %*% t(Lam) + Psi))^2)/length(diag(df.sub$Ge)))
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$`Acc-totGe` <- sqrt(sum((df.sub$Ge[upper.tri(df.sub$Ge,diag = T)] - (Lam %*% t(Lam) + Psi)[upper.tri(df.sub$Ge,diag = T)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$`Acc-totCe` <- sqrt(sum((df.sub$Ce[upper.tri(df.sub$Ge,diag = T)] - cov2cor(Lam %*% t(Lam) + Psi)[upper.tri(df.sub$Ge,diag = T)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$`Acc-0` <- cor(df.sub$trueG,rowMeans(estInt.reg)) # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$`Acc-1` <- cor(df$trueG,rowMeans(estInt.reg)) # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$`Acc-2` <- cor(df$trueG,rowMeans(estInt.reg + estInt.diag))
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$`Acc-3` <- mean(diag(cor(df.sub$trueGE,estInt.reg)))
df.output[df.output$Rep == Rep & df.output$Model == "FA1",]$`Acc-4` <- mean(diag(cor(df.sub$trueGE,estInt.reg + estInt.diag)))

# 5. Factor analytic 2 -----------------
cat("--> FA2 model \n")
# Assign start values
vars = summary(asr.fa1)$varcomp
sv = asreml(y ~ 1 + env,
            random = ~ rr(env, 2):id + 
              diag(env):id,
            residual = ~dsum(~units|env),
            na.action = na.method("include"),
            start.values = TRUE,
            data = df.sub$met.df,
            workspace = "1Gb")
sv = sv$vparameters.table
rownames(sv) <- sv$Component

# Match terms
# rownames(vars)[!rownames(vars) %in% sv$Component]
# sv$Component[!sv$Component %in% rownames(vars)]
# Assign start values
sv[rownames(vars)[rownames(vars)%in%sv$Component],2] <- vars[rownames(vars)%in%sv$Component,1]
sv[grep("env.*id.*env",sv$Component),2] <- 0.6*vars[grep("env.*id.*env",rownames(vars)),1]
sv[grep("rr.*fa1",sv$Component),2] <- vars[grep("rr.*fa1",rownames(vars)),1]

# Model
asr.fa2 <- asreml(y ~ 1 + env,
                  random = ~ rr(env, 2):id + 
                    diag(env):id,
                  residual = ~dsum(~units|env),
                  na.action = na.method("include"),
                  G.param = sv, R.param = sv,
                  data = df.sub$met.df,
                  workspace = "1Gb")
while (asr.fa2$converge != TRUE) {
  asr.fa2 = update.asreml(asr.fa2)
}
end_time <- Sys.time()

# Get effects
estInt.reg  <- matrix(asr.fa2$coefficients$random[grep("^rr.*id", rownames(asr.fa2$coefficients$random))][1:(n_genos*n_envs)],ncol = n_envs,byrow = F)
estInt.diag <- matrix(asr.fa2$coefficients$random[grep("^env.*id", rownames(asr.fa2$coefficients$random))][1:(n_genos*n_envs)],ncol = n_envs,byrow = F)
Lam <- matrix(asr.fa2$vparameters[grep("rr.*fa", names(asr.fa2$vparameters))], ncol = 2)
Psi <- diag(asr.fa2$vparameters[grep("^env.*id", names(asr.fa2$vparameters))])

# Report
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$LogLik  <- asr.fa2$loglik
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$AIC     <- summary.asreml(asr.fa2)$aic
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$BIC     <- summary.asreml(asr.fa2)$bic
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$`runTime` <- difftime(end_time, start_time, units = "secs")
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$VAF     <- sum(diag(Lam %*% t(Lam)))/sum(diag(Lam %*% t(Lam) + Psi))
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$`Acc-varG`  <- sqrt(sum((diag(df.sub$Ge) - diag(Lam %*% t(Lam) + Psi))^2)/length(diag(df.sub$Ge)))
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$`Acc-totGe` <- sqrt(sum((df.sub$Ge[upper.tri(df.sub$Ge,diag = T)] - (Lam %*% t(Lam) + Psi)[upper.tri(df.sub$Ge,diag = T)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$`Acc-totCe` <- sqrt(sum((df.sub$Ce[upper.tri(df.sub$Ge,diag = T)] - cov2cor(Lam %*% t(Lam) + Psi)[upper.tri(df.sub$Ge,diag = T)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$`Acc-0` <- cor(df.sub$trueG,rowMeans(estInt.reg)) # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$`Acc-1` <- cor(df$trueG,rowMeans(estInt.reg)) # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$`Acc-2` <- cor(df$trueG,rowMeans(estInt.reg + estInt.diag))
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$`Acc-3` <- mean(diag(cor(df.sub$trueGE,estInt.reg)))
df.output[df.output$Rep == Rep & df.output$Model == "FA2",]$`Acc-4` <- mean(diag(cor(df.sub$trueGE,estInt.reg + estInt.diag)))

# 6. Factor analytic 3 -----------------
cat("--> FA3 model \n")
# Assign start values
vars = summary(asr.fa2)$varcomp
sv = asreml(y ~ 1 + env,
            random = ~ rr(env, 3):id + 
              diag(env):id,
            residual = ~dsum(~units|env),
            na.action = na.method("include"),
            start.values = TRUE,
            data = df.sub$met.df,
            workspace = "1Gb")
sv = sv$vparameters.table
rownames(sv) <- sv$Component

# Assign start values
sv[rownames(vars)[rownames(vars)%in%sv$Component],2] <- vars[rownames(vars)%in%sv$Component,1]
sv[grep("env.*id.*env",sv$Component),2] <- 0.6*vars[grep("env.*id.*env",rownames(vars)),1]
sv[grep("rr.*fa1",sv$Component),2] <- vars[grep("rr.*fa1",rownames(vars)),1]
sv[grep("rr.*fa2",sv$Component),2] <- vars[grep("rr.*fa2",rownames(vars)),1]

# Model
asr.fa3 <- asreml(y ~ 1 + env,
                  random = ~ rr(env, 3):id + 
                    diag(env):id,
                  residual = ~dsum(~units|env),
                  na.action = na.method("include"),
                  G.param = sv, R.param = sv,
                  data = df.sub$met.df,
                  workspace = "1Gb")
while (asr.fa3$converge != TRUE) {
  asr.fa3 = update.asreml(asr.fa3)
}
end_time <- Sys.time()

# Get effects
estInt.reg  <- matrix(asr.fa3$coefficients$random[grep("^rr.*id", rownames(asr.fa3$coefficients$random))][1:(n_genos*n_envs)],ncol = n_envs,byrow = F)
estInt.diag <- matrix(asr.fa3$coefficients$random[grep("^env.*id", rownames(asr.fa3$coefficients$random))][1:(n_genos*n_envs)],ncol = n_envs,byrow = F)
Lam <- matrix(asr.fa3$vparameters[grep("rr.*fa", names(asr.fa3$vparameters))], ncol = 3)
Psi <- diag(asr.fa3$vparameters[grep("^env.*id", names(asr.fa3$vparameters))])
# Report
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$LogLik  <- asr.fa3$loglik
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$AIC     <- summary.asreml(asr.fa3)$aic
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$BIC     <- summary.asreml(asr.fa3)$bic
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$`runTime` <- difftime(end_time, start_time, units = "secs")
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$VAF     <- sum(diag(Lam %*% t(Lam)))/sum(diag(Lam %*% t(Lam) + Psi))
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$`Acc-varG`  <- sqrt(sum((diag(df.sub$Ge) - diag(Lam %*% t(Lam) + Psi))^2)/length(diag(df.sub$Ge)))
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$`Acc-totGe` <- sqrt(sum((df.sub$Ge[upper.tri(df.sub$Ge,diag = T)] - (Lam %*% t(Lam) + Psi)[upper.tri(df.sub$Ge,diag = T)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$`Acc-totCe` <- sqrt(sum((df.sub$Ce[upper.tri(df.sub$Ge,diag = T)] - cov2cor(Lam %*% t(Lam) + Psi)[upper.tri(df.sub$Ge,diag = T)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$`Acc-0` <- cor(df.sub$trueG,rowMeans(estInt.reg)) # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$`Acc-1` <- cor(df$trueG,rowMeans(estInt.reg)) # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$`Acc-2` <- cor(df$trueG,rowMeans(estInt.reg + estInt.diag))
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$`Acc-3` <- mean(diag(cor(df.sub$trueGE,estInt.reg)))
df.output[df.output$Rep == Rep & df.output$Model == "FA3",]$`Acc-4` <- mean(diag(cor(df.sub$trueGE,estInt.reg + estInt.diag)))

# 7. Factor analytic 4 -----------------
cat("--> FA4 model \n")
# Assign start values
vars = summary(asr.fa3)$varcomp
sv = asreml(y ~ 1 + env,
            random = ~ rr(env, 4):id + 
              diag(env):id,
            residual = ~dsum(~units|env),
            na.action = na.method("include"),
            start.values = TRUE,
            data = df.sub$met.df,
            workspace = "1Gb")
sv = sv$vparameters.table
rownames(sv) <- sv$Component

# Assign start values
sv[rownames(vars)[rownames(vars)%in%sv$Component],2] <- vars[rownames(vars)%in%sv$Component,1]
sv[grep("env.*id.*env",sv$Component),2] <- 0.6*vars[grep("env.*id.*env",rownames(vars)),1]
sv[grep("rr.*fa1",sv$Component),2] <- vars[grep("rr.*fa1",rownames(vars)),1]
sv[grep("rr.*fa2",sv$Component),2] <- vars[grep("rr.*fa2",rownames(vars)),1]
sv[grep("rr.*fa3",sv$Component),2] <- vars[grep("rr.*fa3",rownames(vars)),1]

# Model
asr.fa4 <- asreml(y ~ 1 + env,
                  random = ~ rr(env, 4):id + 
                    diag(env):id,
                  residual = ~dsum(~units|env),
                  na.action = na.method("include"),
                  G.param = sv, R.param = sv,
                  data = df.sub$met.df,
                  workspace = "1Gb")
while (asr.fa4$converge != TRUE) {
  asr.fa4 = update.asreml(asr.fa4)
}
# Reinsert fa4 values
vars = summary(asr.fa4)$varcomp
sv[rownames(vars)[rownames(vars)%in%sv$Component],2] <- vars[rownames(vars)%in%sv$Component,1]
asr.fa4 <- asreml(y ~ 1 + env,
                  random = ~ rr(env, 4):id + 
                    diag(env):id,
                  residual = ~dsum(~units|env),
                  na.action = na.method("include"),
                  G.param = sv, R.param = sv,
                  data = df.sub$met.df,
                  workspace = "1Gb")
while (asr.fa4$converge != TRUE) {
  asr.fa4 = update.asreml(asr.fa4)
}
end_time <- Sys.time()

# Get effects
estInt.reg  <- matrix(asr.fa4$coefficients$random[grep("^rr.*id", rownames(asr.fa4$coefficients$random))][1:(n_genos*n_envs)],ncol = n_envs,byrow = F)
estInt.diag <- matrix(asr.fa4$coefficients$random[grep("^env.*id", rownames(asr.fa4$coefficients$random))][1:(n_genos*n_envs)],ncol = n_envs,byrow = F)
Lam <- matrix(asr.fa4$vparameters[grep("rr.*fa", names(asr.fa4$vparameters))], ncol = 4)
Psi <- diag(asr.fa4$vparameters[grep("^env.*id", names(asr.fa4$vparameters))])
# Report
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$LogLik  <- asr.fa4$loglik
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$AIC     <- summary.asreml(asr.fa4)$aic
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$BIC     <- summary.asreml(asr.fa4)$bic
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$`runTime` <- difftime(end_time, start_time, units = "secs")
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$VAF     <- sum(diag(Lam %*% t(Lam)))/sum(diag(Lam %*% t(Lam) + Psi))
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$`Acc-varG`  <- sqrt(sum((diag(df.sub$Ge) - diag(Lam %*% t(Lam) + Psi))^2)/length(diag(df.sub$Ge)))
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$`Acc-totGe` <- sqrt(sum((df.sub$Ge[upper.tri(df.sub$Ge,diag = T)] - (Lam %*% t(Lam) + Psi)[upper.tri(df.sub$Ge,diag = T)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$`Acc-totCe` <- sqrt(sum((df.sub$Ce[upper.tri(df.sub$Ge,diag = T)] - cov2cor(Lam %*% t(Lam) + Psi)[upper.tri(df.sub$Ge,diag = T)])^2)/length(df.sub$Ge[upper.tri(df.sub$Ge,F)]))
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$`Acc-0` <- cor(df.sub$trueG,rowMeans(estInt.reg)) # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$`Acc-1` <- cor(df$trueG,rowMeans(estInt.reg)) # implicit main effect
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$`Acc-2` <- cor(df$trueG,rowMeans(estInt.reg + estInt.diag))
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$`Acc-3` <- mean(diag(cor(df.sub$trueGE,estInt.reg)))
df.output[df.output$Rep == Rep & df.output$Model == "FA4",]$`Acc-4` <- mean(diag(cor(df.sub$trueGE,estInt.reg + estInt.diag)))

# Check results -----------------
df.output

# Table contains columns for:
# Column       Description
# Scenario:    Name of scenario
# subScenario: Possible subscenario (in this case a scenario without spatial variation was run)
# Envs:        Size of the MET dataset in terms of environments.
# Rep:         Replicate number
# Model:       Type of statistical model
# LogLik:      Log likelihood 
# AIC:         Akaike Infomation Criteria
# BIC:         Bayesian Information Criteria
# runTime:     Time for model to converge
# VAF:         Proportion of genetic variance explained by the model
# Acc-varG:    Correlation between true and estimated environment genetic variances
# Acc-totGe:   Correlation between true and estimated between-environment genetic variance matrix
# Acc-totCe:   Correlation between true and estimated between-environment genetic correlation matrix
# Acc-0:       Correlation between true and estimated main (additive) genetic effects in MET
# Acc-1:       Correlation between true additive main genetic effects in TPE and estimated by model in MET
# Acc-2:       Correlation between true total main genetic effects in TPE and estimated by model in MET
# Acc-3:       Correlation between true additive genotype by environment effects in TPE and estimated by model in MET
# Acc-4:       Correlation between true total genotype by environment effects in TPE and estimated by model in MET
# varG:        Proportion of variance explained by genetic main effect.
# varGEI:      Proportion of variance explained by genotype by environment interaction.
# nGEI:        Proportion of variance explained by non-crossover genotype by environment interaction.
# cGEI:        Proportion of variance explained by crossover genotype by environment interaction.

# Save results -----------------
cat("Saving rep:", Rep, "\n")
file.name <- paste0("Results_",scenario,subscenario,n_envs,".csv")
write.table(df.output, file.name, sep = ",", col.names = !file.exists(file.name), row.names = F, append = T)
