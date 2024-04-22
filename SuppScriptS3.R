
##########################################################
#
# Supplementary Script S3: A framework for simulating GEI
# 
# Comparison of models fitted to a simulated MET dataset
#
# Script authors: J. Bancic and D.J. Tolhurst
#
##########################################################

# This script briefly demonstrates model comparison using ASReml-R

# It is important to note that the models here are intended 
# for demonstration purposes only, and should be tailored to
# the research objectives. 

# Load ASReml-R
library(asreml)

# The following are generated in Supplementary Script S2
met_df   # MET data frame 
gvs_met  # true genetic values in MET
gvs_tpe  # true genetic values in MET

# Output data frame
modelNames <- c("Maximum", "Expected", "Main", "Comp", "MDiag", "Diag", "FA1")
output <- data.frame(
  GEI = rep("Moderate", length(modelNames)),
  Envs = rep(pm, length(modelNames)),
  Model = modelNames,
  LogLik = NA,
  AIC = NA,
  acc_g_met = NA,
  acc_g_tpe = NA,
  acc_ge_met = NA,
  stringsAsFactors = FALSE
)
output

# Get true values
g_met_true <- with(gvs_met, tapply(gv.Trait1, id, mean))
g_tpe_true <- with(gvs_tpe, tapply(gv.Trait1, id, mean))
ge_met_true <- with(gvs_met, tapply(gv.Trait1, list(id, env), mean))
(Ge_vars <- measure_variances(Ge))

#-- Expected accuracies
# MET-TPE alignment
output[output$Model == "Maximum",]$acc_g_tpe <- 
  sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/pm))
# Main effect accuracy in TPE
output[output$Model == "Expected",]$acc_g_tpe <- 
  sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/pm + sigm2e/pm/r))
# Main effect accuracy in MET
output[output$Model == "Expected",]$acc_g_met <- 
  sqrt((Ge_vars[1,2] + Ge_vars[2,2]/pm)/(Ge_vars[1,2] + Ge_vars[2,2]/pm + sigm2e/pm/r))
# Accuracy for GE effects in MET
output[output$Model == "Expected",]$acc_ge_met <- 
  sqrt((Ge_vars[1,2] + Ge_vars[2,2])/(Ge_vars[1,2] + Ge_vars[2,2] + sigm2e/r))

#-- Model 1: Main effects only
asr.main <- asreml(y.Trait1 ~ 1 + env,
                   random   = ~ id + diag(env):block,
                   residual = ~ dsum(~ ar1(col):ar1(row) | env),
                   data     = met_df, 
                   workspace = "1Gb")
while (asr.main$converge != TRUE) {asr.main = update.asreml(asr.main)}

# Predicted genotype main effects
g_met_pred = asr.main$coefficients$random[
  grep(pattern = "^id", rownames(asr.main$coefficients$random)),]

# Store output
cat("LogLik: ", output[output$Model == "Main",]$LogLik <- asr.main$loglik,"\n")
cat("AIC: ",    output[output$Model == "Main",]$AIC <- summary.asreml(asr.main)$aic,"\n")
cat("acc_g_met: ", output[output$Model == "Main",]$acc_g_met <- cor(g_met_pred,g_met_true),"\n")
cat("acc_g_tpe: ", output[output$Model == "Main",]$acc_g_tpe <- cor(g_met_pred,g_tpe_true),"\n")
cat("acc_ge_met: ", output[output$Model == "Main",]$acc_ge_met <- NA)


#-- Model 2: Compound symmetry
asr.comp <- asreml(y.Trait1 ~ 1 + env,
                   random   = ~ id + env:id + diag(env):block,
                   residual = ~ dsum(~ ar1(col):ar1(row) | env),
                   data     = met_df,
                   workspace = "1Gb")
while (asr.comp$converge != TRUE) {asr.comp = update.asreml(asr.comp)}

# Predicted genotype main effects
g_met_pred = asr.comp$coefficients$random[
  grep(pattern = "^id", rownames(asr.comp$coefficients$random)),]
ge_met_pred = asr.comp$coefficients$random[
  grep(pattern = "^env.*id", rownames(asr.comp$coefficients$random)),]
ge_met_pred = matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)

# Store output
cat("LogLik: ", output[output$Model == "Comp",]$LogLik <- asr.comp$loglik,"\n")
cat("AIC: ",    output[output$Model == "Comp",]$AIC <- summary.asreml(asr.comp)$aic,"\n")
cat("acc_g_met: ", output[output$Model == "Comp",]$acc_g_met <- cor(g_met_pred,g_met_true),"\n")
cat("acc_g_tpe: ", output[output$Model == "Comp",]$acc_g_tpe <- cor(g_met_pred,g_tpe_true),"\n")
cat("acc_ge_met: ", output[output$Model == "Comp",]$acc_ge_met <- mean(diag(cor(ge_met_pred,ge_met_true))),"\n")


#-- Model 3: Main effects + diagonal
asr.mdiag <- asreml(y.Trait1 ~ 1 + env,
                    random   = ~ id + diag(env):id + diag(env):block,
                    residual = ~ dsum(~ ar1(col):ar1(row) | env),
                    data     = met_df,
                    workspace = "1Gb")
while (asr.mdiag$converge != TRUE) {asr.mdiag = update.asreml(asr.mdiag)}

# Predicted genotype main effects
g_met_pred = asr.mdiag$coefficients$random[
  grep(pattern = "^id", rownames(asr.mdiag$coefficients$random)),]
ge_met_pred = asr.mdiag$coefficients$random[
  grep(pattern = "^env.*id", rownames(asr.mdiag$coefficients$random)),]
ge_met_pred = matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)

# Store output
cat("LogLik: ", output[output$Model == "MDiag",]$LogLik <- asr.mdiag$loglik,"\n")
cat("AIC: ",    output[output$Model == "MDiag",]$AIC <- summary.asreml(asr.mdiag)$aic,"\n")
cat("acc_g_met: ", output[output$Model == "MDiag",]$acc_g_met <- cor(g_met_pred,g_met_true),"\n")
cat("acc_g_tpe: ", output[output$Model == "MDiag",]$acc_g_tpe <- cor(g_met_pred,g_tpe_true),"\n")
cat("acc_ge_met: ", output[output$Model == "MDiag",]$acc_ge_met <- mean(diag(cor(ge_met_pred,ge_met_true))),"\n")


#-- Model 4: Diagonal
asr.diag <- asreml(y.Trait1 ~ 1 + env,
                   random   = ~ diag(env):id + diag(env):block,
                   residual = ~ dsum(~ ar1(col):ar1(row) | env),
                   data = met_df,
                   workspace = "1Gb")
while (asr.diag$converge != TRUE) {asr.diag = update.asreml(asr.diag)}

# Predicted genotype main effects
ge_met_pred = asr.diag$coefficients$random[
  grep(pattern = "^env.*id", rownames(asr.diag$coefficients$random)),]
ge_met_pred = matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)

# Store output
cat("LogLik: ", output[output$Model == "Diag",]$LogLik <- asr.diag$loglik,"\n")
cat("AIC: ",    output[output$Model == "Diag",]$AIC <- summary.asreml(asr.diag)$aic,"\n")
cat("acc_g_met: ", output[output$Model == "Diag",]$acc_g_met <- NA,"\n")
cat("acc_g_tpe: ", output[output$Model == "Diag",]$acc_g_tpe <- NA,"\n")
cat("acc_ge_met: ", output[output$Model == "Diag",]$acc_ge_met <- mean(diag(cor(ge_met_pred,ge_met_true))),"\n")


#-- Model 5: Factor analytic of order 1
asr.fa1 <- asreml(y.Trait1 ~ 1 + env,
                  random   = ~ rr(env, 1):id + diag(env):id + diag(env):block,
                  residual = ~ dsum(~ ar1(col):ar1(row) | env),
                  data = met_df,
                  workspace = "2Gb")
while (asr.fa1$converge != TRUE) {asr.fa1 = update.asreml(asr.fa1)}

# Predicted genotype main effects
ge_met_pred = asr.fa1$coefficients$random[
  grep("^rr.*id", rownames(asr.fa1$coefficients$random)),]
ge_met_pred = ge_met_pred[
  grep("Comp", names(ge_met_pred), invert = TRUE)]
ge_met_pred = matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)
g_met_pred = rowMeans(ge_met_pred)

# Store output
cat("LogLik: ", output[output$Model == "FA1",]$LogLik <- asr.fa1$loglik,"\n")
cat("AIC: ",    output[output$Model == "FA1",]$AIC <- summary.asreml(asr.fa1)$aic,"\n")
cat("acc_g_met: ", output[output$Model == "FA1",]$acc_g_met <- cor(g_met_pred,g_met_true),"\n")
cat("acc_g_tpe: ", output[output$Model == "FA1",]$acc_g_tpe <- cor(g_met_pred,g_tpe_true),"\n")
cat("acc_ge_met: ", output[output$Model == "FA1",]$acc_ge_met <- mean(diag(cor(ge_met_pred,ge_met_true))),"\n")

# View output example
output
#        GEI Envs    Model    LogLik      AIC acc_g_met acc_g_tpe acc_ge_met
# 1 Moderate   20  Maximum        NA       NA        NA 0.9501404         NA
# 2 Moderate   20 Expected        NA       NA 0.9158263 0.8701636  0.6522638
# 3 Moderate   20     Main -19457.15 39076.29 0.9348395 0.8881698         NA
# 4 Moderate   20     Comp -18910.10 37984.21 0.9382199 0.8729845  0.5631765
# 5 Moderate   20    MDiag -18653.99 37509.98 0.9275970 0.8909391  0.5726093
# 6 Moderate   20     Diag -18924.11 38048.22        NA        NA  0.6924883
# 7 Moderate   20      FA1 -18534.16 37308.33 0.8942324 0.8798488  0.5369332

# The above can be extended to higher order factor analytic models by changing rr(env, 1):id accordingly


# end of script


