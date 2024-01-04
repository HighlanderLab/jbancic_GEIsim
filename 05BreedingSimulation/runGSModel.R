########################################################################
## Run genomic model
########################################################################
cat(" Running GS model\n")
# NOTE: Use records from PYT, AYT and EYT trials to predict performance 
# of DHs, and perform genomic assisted selection of AYTs and EYTs

#-----------------------------------------------------------------------
# Prepare data to run the model
#-----------------------------------------------------------------------
# Add genotypes with missing phenotypes
tmpMET = MET # specify training population

# Square: Subset to get only the last 1 year of data
tmpMET = tmpMET[tmpMET$Year >= (year-1),]

# Triangle: Including complete contemporary groups
# tmpMET <- tmpMET[!(tmpMET$Stage == "EYT" & tmpMET$Year == (year-2)),]
# tmpMET <- tmpMET[!(tmpMET$Stage == "AYT" & tmpMET$Year == (year-3)),]
# tmpMET <- tmpMET[!(tmpMET$Stage == "EYT" & tmpMET$Year == (year-3)),]

# Turn variables into factors in MET training dataset
tmpMET$Stage  = factor(tmpMET$Stage,  levels = unique(as.character(tmpMET$Stage)))
tmpMET$Loc    = factor(tmpMET$Loc,    levels = unique(as.character(tmpMET$Loc)))
tmpMET$Year   = factor(tmpMET$Year,   levels = unique(as.character(tmpMET$Year)))
tmpMET$Env    = factor(paste0(tmpMET$Year,"-",tmpMET$Loc))
tmpMET$Id     = factor(tmpMET$Id,     levels = unique(as.character(tmpMET$Id)))
tmpMET$Male   = factor(tmpMET$Male,   levels = unique(as.character(tmpMET$Male)))
tmpMET$Female = factor(tmpMET$Female, levels = unique(as.character(tmpMET$Female)))
tmpMET$Comb   = factor(tmpMET$Comb,   levels = unique(as.character(tmpMET$Comb)))
tmpMET        = tmpMET[order(tmpMET$Env, tmpMET$Loc),]
cat(" -> MET includes",length(levels(tmpMET$Year)), "year(s) with", dim(tmpMET)[1], 
    "observations,", length(levels(tmpMET$Id)), "genotypes,", "and stages", levels(tmpMET$Stage),"\n")

# Get marker matrix of tested genotypes
# Mtt = pullSnpGeno(GenoPop[which(GenoPop@id %in% tmpMET$Id == T)])-1 # specify pop
Mtt = GenoPop[which(rownames(GenoPop) %in% tmpMET$Id == T),]-1
Mtt = Mtt[levels(tmpMET$Id),]

# Create GRM with simple scaling (to help with convergence)
Ktt = (Mtt %*% t(Mtt))/10000
diag(Ktt) = diag(Ktt) + 0.00001

# Get an inverse of K_tt and input format
Ktti = G.inverse(G = Ktt, sparseform = T)$Ginv

#-----------------------------------------------------------------------
# Run main effects plus diagonal model
#-----------------------------------------------------------------------
if (!exists("GSmodel")) {
  # Run model the first year
  GSmodel = asreml(fixed     = Pheno ~ 1 + Env,
                   random    = ~ vm(Id, Ktti) + idv(Env):vm(Id, Ktti),
                   # diag(Stage):Env, # not necessary
                   residual  = ~ units,
                   na.action = na.method("include"),
                   data      = tmpMET,
                   workspace = "5e8")
  # Rerun if not converged
  while (GSmodel$converge != TRUE) {
    GSmodel = update.asreml(GSmodel)
  }
} else {
  # Take variance parameters from previous year model as start values
  vars = summary(GSmodel)$varcomp
  GSmodel.sv = asreml(fixed     = Pheno ~ 1 + Env,
                      random    = ~ vm(Id, Ktti) + idv(Env):vm(Id, Ktti),
                      residual  = ~ units,
                      na.action = na.method("include"),
                      start.values = TRUE,
                      data      = tmpMET,
                      workspace = "5e8")
  GSmodel.sv = GSmodel.sv$vparameters.table
  # Assign start values
  take.sv <- which(!is.na(match(GSmodel.sv$Component,rownames(vars))))
  take.vars <- match(GSmodel.sv$Component,rownames(vars))[
    !is.na(match(GSmodel.sv$Component,rownames(vars)))]
  GSmodel.sv[take.sv,2] <- vars[take.vars,1]
  
  # Run model
  GSmodel = asreml(fixed     = Pheno ~ 1 + Env,
                   random    = ~ vm(Id, Ktti) + idv(Env):vm(Id, Ktti),
                   # diag(Stage):Env, # not necessary
                   residual  = ~ units,
                   na.action = na.method("include"),
                   G.param   = GSmodel.sv,
                   R.param   = GSmodel.sv, 
                   data      = tmpMET,
                   workspace = "5e8")
  # Rerun if not converged
  while (GSmodel$converge != TRUE) {
    GSmodel = update.asreml(GSmodel)
  }
}

## Get main genotype effects
temp  = GSmodel$coefficients$random[
  grep(pattern = ".*Id", rownames(GSmodel$coefficients$random)),]
MainG = temp[grep(pattern = "Env", names(temp),invert = T)]
names(MainG) = sub(pattern = ".*_","", names(MainG))
# Add intercept back in
MainG = MainG + GSmodel$coefficients$fixed[
  grep(pattern = "Intercept", rownames(GSmodel$coefficients$fixed)),]

#-----------------------------------------------------------------------
# Assign estimated additive effects for selection
#-----------------------------------------------------------------------
cat(" Predicting HDRWs \n")
# NOTE: Make sure you know where the values are assigned to!
# Predict additive effects of untested HDRWs
Mut   = pullSnpGeno(HDRW)-1 # specify untested pop
G12   = Mtt %*% t(Mut)/10000; dim(G12) 
ypred = t(G12) %*% solve(Ktt) %*% MainG
cat("  DH-MET-PS:", cor(getPopParams(HDRW)$Mean_MET, getPopParams(HDRW)$estMean))
cat("  DH-TPE-PS:", cor(getPopParams(HDRW)$Mean_TPE, getPopParams(HDRW)$estMean), "\n")
HDRW  = setMisc(HDRW, "estMean", ypred[which(HDRW@id %in% rownames(ypred))]) 
cat("  DH-MET-GS:", cor(getPopParams(HDRW)$Mean_MET, getPopParams(HDRW)$estMean))
cat("  DH-TPE-GS:", cor(getPopParams(HDRW)$Mean_TPE, getPopParams(HDRW)$estMean), "\n")

cat(" Performing genomic-assisted selection \n")
# Set genomic-estimated additive effects for PYTs
temp = MainG[which(names(MainG) %in% PYT@id == T)]
cat("  PYT-MET-PS:", cor(getPopParams(PYT)$Mean_MET,getPopParams(PYT)$estMean))
cat("  PYT-TPE-PS:", cor(getPopParams(PYT)$Mean_TPE,getPopParams(PYT)$estMean), "\n")
PYT  = setMisc(PYT, "estMean", temp[PYT@id])
cat("  PYT-MET-GS:", cor(getPopParams(PYT)$Mean_MET,getPopParams(PYT)$estMean))
cat("  PYT-TPE-GS:", cor(getPopParams(PYT)$Mean_TPE,getPopParams(PYT)$estMean), "\n")

# Set genomic-estimated additive effects for AYTs
temp = MainG[which(names(MainG) %in% AYT@id == T)]
cat("  AYT-MET-PS:", cor(getPopParams(AYT)$Mean_MET,getPopParams(AYT)$estMean), "\n") 
AYT  = setMisc(AYT, "estMean", temp[AYT@id])
cat("  AYT-MET-GS:", cor(getPopParams(AYT)$Mean_MET,getPopParams(AYT)$estMean), "\n") 

#-----------------------------------------------------------------------
## Clean up the environment
rm(tmpMET, Mtt, Ktt, Ktti, temp, MainG, Mut, G12)
gc(verbose = F)