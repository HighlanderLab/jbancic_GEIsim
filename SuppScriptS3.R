
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
pe <- length(levels(gvs_met$env))

# Output data frame
modelNames <- c("MET-TPE", "Expected", "Main", "Comp", "MDiag", "Diag", "FA1")
output <- data.frame(
  GEI = rep("Moderate", length(modelNames)),
  Envs = rep(pm, length(modelNames)),
  Model = modelNames,
  LogLik = NA,
  AIC  = NA,
  r_g  = NA,
  r_m  = NA,
  r_ge = NA,
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
output[output$Model == "MET-TPE",]$r_g <- 
  sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/pm))
# Main effect accuracy in TPE
output[output$Model == "Expected",]$r_g <- 
  sqrt(Ge_vars[1,2]/(Ge_vars[1,2] + Ge_vars[2,2]/pm + sigm2e/pm/r))
# Main effect accuracy in MET
output[output$Model == "Expected",]$r_m <- 
  sqrt((Ge_vars[1,2] + Ge_vars[2,2]/pm)/(Ge_vars[1,2] + Ge_vars[2,2]/pm + sigm2e/pm/r))
# Accuracy for GE effects in MET
output[output$Model == "Expected",]$r_ge <- 
  sqrt((Ge_vars[1,2] + Ge_vars[2,2])/(Ge_vars[1,2] + Ge_vars[2,2] + sigm2e/r))



#-- Model 1: Main effects only
asr.main <- asreml(y.Trait1 ~ 1 + env,
                   random   = ~ id + diag(env):block,
                   residual = ~ dsum(~ ar1(col):ar1(row) | env),
                   data     = met_df, 
                   workspace = "1Gb")
while (asr.main$converge != TRUE) {asr.main = update.asreml(asr.main)}

# Predicted genotype main effects
g_met_pred <- asr.main$coefficients$random[
  grep(pattern = "^id", rownames(asr.main$coefficients$random)),]
# Repeated main effects across all environments to obtain genotype by environment effects
ge_met_pred <- matrix(g_met_pred, ncol = pe, nrow = v)
# Total genetic effects not calculated

# Store output
cat("LogLik: ", output[output$Model == "Main",]$LogLik <- asr.main$loglik,"\n")
cat("AIC: ",    output[output$Model == "Main",]$AIC <- summary.asreml(asr.main)$aic,"\n")
cat("r_g: ",  output[output$Model == "Main",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_m: ",  output[output$Model == "Main",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_ge: ", output[output$Model == "Main",]$r_ge <- mean(diag(cor(ge_met_pred,ge_met_true))),"\n")



#-- Model 2: Compound symmetry
asr.comp <- asreml(y.Trait1 ~ 1 + env,
                   random   = ~ id + env:id + diag(env):block,
                   residual = ~ dsum(~ ar1(col):ar1(row) | env),
                   data     = met_df,
                   workspace = "1Gb")
while (asr.comp$converge != TRUE) {asr.comp = update.asreml(asr.comp)}

# Predicted genotype main effects
g_met_pred <- asr.comp$coefficients$random[
  grep(pattern = "^id", rownames(asr.comp$coefficients$random)),]
# Predicted genotype by environment effects
ge_met_pred <- asr.comp$coefficients$random[
  grep(pattern = "^env.*id", rownames(asr.comp$coefficients$random)),]
ge_met_pred <- matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)
# Predicted total genotype by environment effects
ge_tot_met_pred <- g_met_pred + ge_met_pred

# Store output
cat("LogLik: ", output[output$Model == "Comp",]$LogLik <- asr.comp$loglik,"\n")
cat("AIC: ",    output[output$Model == "Comp",]$AIC <- summary.asreml(asr.comp)$aic,"\n")
cat("r_m: ",  output[output$Model == "Comp",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_g: ",  output[output$Model == "Comp",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_ge: ", output[output$Model == "Comp",]$r_ge <- mean(diag(cor(ge_tot_met_pred,ge_met_true))),"\n")



#-- Model 3: Main effects + diagonal
asr.mdiag <- asreml(y.Trait1 ~ 1 + env,
                    random   = ~ id + diag(env):id + diag(env):block,
                    residual = ~ dsum(~ ar1(col):ar1(row) | env),
                    data     = met_df,
                    workspace = "1Gb")
while (asr.mdiag$converge != TRUE) {asr.mdiag = update.asreml(asr.mdiag)}

# Predicted genotype main effects
g_met_pred <- asr.mdiag$coefficients$random[
  grep(pattern = "^id", rownames(asr.mdiag$coefficients$random)),]
# Predicted genotype by environment effects
ge_met_pred <- asr.mdiag$coefficients$random[
  grep(pattern = "^env.*id", rownames(asr.mdiag$coefficients$random)),]
ge_met_pred <- matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)
# Predicted total genotype by environment effects
ge_tot_met_pred <- g_met_pred + ge_met_pred

# Store output
cat("LogLik: ", output[output$Model == "MDiag",]$LogLik <- asr.mdiag$loglik,"\n")
cat("AIC: ",    output[output$Model == "MDiag",]$AIC <- summary.asreml(asr.mdiag)$aic,"\n")
cat("r_m: ",  output[output$Model == "MDiag",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_g: ",  output[output$Model == "MDiag",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_ge: ", output[output$Model == "MDiag",]$r_ge <- mean(diag(cor(ge_tot_met_pred,ge_met_true))),"\n")



#-- Model 4: Diagonal
asr.diag <- asreml(y.Trait1 ~ 1 + env,
                   random   = ~ diag(env):id + diag(env):block,
                   residual = ~ dsum(~ ar1(col):ar1(row) | env),
                   data = met_df,
                   workspace = "1Gb")
while (asr.diag$converge != TRUE) {asr.diag = update.asreml(asr.diag)}

# Predicted genotype by environment effects
ge_met_pred <- asr.diag$coefficients$random[
  grep(pattern = "^env.*id", rownames(asr.diag$coefficients$random)),]
ge_met_pred <- matrix(ge_met_pred, nrow = v, ncol = pm, byrow = F)
# Implicit main genotype effects
g_met_pred <- rowMeans(ge_met_pred)
# Total genetic effects not calculated

# Store output
cat("LogLik: ", output[output$Model == "Diag",]$LogLik <- asr.diag$loglik,"\n")
cat("AIC: ",    output[output$Model == "Diag",]$AIC <- summary.asreml(asr.diag)$aic,"\n")
cat("r_g: ",  output[output$Model == "Diag",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_m: ",  output[output$Model == "Diag",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_ge: ", output[output$Model == "Diag",]$r_ge <- mean(diag(cor(ge_met_pred,ge_met_true))),"\n")



#-- Model 5: Factor analytic of order 1
asr.fa1 <- asreml(y.Trait1 ~ 1 + env,
                  random   = ~ rr(env, 1):id + diag(env):id + diag(env):block,
                  residual = ~ dsum(~ ar1(col):ar1(row) | env),
                  data = met_df,
                  workspace = "2Gb")
while (asr.fa1$converge != TRUE) {asr.fa1 = update.asreml(asr.fa1)}

# Predicted regression genotype by environment effects
ge.reg_met_pred <- asr.fa1$coefficients$random[
  grep("^rr.*id", rownames(asr.fa1$coefficients$random)),]
ge.reg_met_pred <- ge.reg_met_pred[
  grep("Comp", names(ge.reg_met_pred), invert = TRUE)]
ge.reg_met_pred <- matrix(ge.reg_met_pred, nrow = v, ncol = pm, byrow = F)
# Predicted residual genotype by environment effects
ge.diag_met_pred <- asr.fa1$coefficients$random[
  grep("^env.*id", rownames(asr.fa1$coefficients$random)),]
ge.diag_met_pred <- matrix(ge.diag_met_pred, nrow = v, ncol = pm, byrow = F)
# Implicit genotype main effects
g_met_pred <- rowMeans(ge.reg_met_pred)
# Predicted total genotype by environment effects
ge_tot_met_pred <- g_met_pred + ge.diag_met_pred

# Store output
cat("LogLik: ", output[output$Model == "FA1",]$LogLik <- asr.fa1$loglik,"\n")
cat("AIC: ",    output[output$Model == "FA1",]$AIC <- summary.asreml(asr.fa1)$aic,"\n")
cat("r_m: ",  output[output$Model == "FA1",]$r_m <- cor(g_met_pred,g_met_true),"\n")
cat("r_g: ",  output[output$Model == "FA1",]$r_g <- cor(g_met_pred,g_tpe_true),"\n")
cat("r_ge: ", output[output$Model == "FA1",]$r_ge <- mean(diag(cor(ge_tot_met_pred,ge_met_true))),"\n")
# Note: The above can be extended to higher order factor analytic models by changing rr(env, 1):id accordingly



#-- View output example
output
#        GEI Envs    Model    LogLik      AIC       r_m       r_g      r_ge
# 1 Moderate   20  MET-TPE        NA       NA 0.9501404        NA        NA
# 2 Moderate   20 Expected        NA       NA 0.8886864 0.9353211 0.7062228
# 3 Moderate   20     Main -16786.62 33735.25 0.9012443 0.9297892 0.5860806
# 4 Moderate   20     Comp -16266.80 32697.61 0.9140460 0.9406453 0.7364472
# 5 Moderate   20    MDiag -16015.77 32233.53 0.9010404 0.9251952 0.7606309
# 6 Moderate   20     Diag -16406.36 33012.72 0.9289324 0.9503625 0.6713395
# 7 Moderate   20      FA1 -15823.94 31887.89 0.8984942 0.9296549 0.7546676

# end of script


