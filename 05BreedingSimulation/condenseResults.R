# Condense results of the simulation by obtaining the means and variances
# rather than storing entire populations

# Make some variables numeric
temp <- trackParams %>% mutate_at(c(1:2), as.numeric)

# Calculate accuracies per stage and melt into a single column
cors <- temp %>%
  group_by(Year, Stage, Scenario, subScenario) %>%
  summarize(acc_Mean_METtoTPE    = cor(Mean_MET, Mean_TPE),    # MET-TPE alignment for GS
            acc_Mean_subMETtoTPE = cor(Mean_subMET, Mean_TPE), # MET-TPE alignment for Pheno
            acc_Mean_subMETtoMET = cor(Mean_subMET, Mean_MET), # true accuracy in MET for Pheno
            acc_Mean_estMETtoTPE = cor(estMean_MET, Mean_TPE), # realised accuaracy in TPE
            acc_Mean_estMETtoMET = cor(estMean_MET, Mean_MET), # realised accuracy in MET
            acc_SD_METtoTPE    = cor(SD_MET, SD_TPE),
            acc_SD_subMETtoTPE = cor(SD_subMET, SD_TPE),
            acc_SD_subMETtoMET = cor(SD_subMET, SD_MET),
            acc_SD_estMETtoMET = cor(estSD_MET, SD_MET),
            acc_SD_setMETtoTPE = cor(estSD_MET, SD_TPE))
cors <- melt(cors,
             id.vars = c("Year","Stage","Scenario","subScenario"),
             measure.vars = colnames(cors)[5:dim(cors)[2]])
colnames(cors)[6] <- "value_mean"

# Melt all other variables into a single column
takeCols <- colnames(temp)[7:length(colnames(temp))]
temp <- melt(temp,
             id.vars = c("Year","Rep","Stage","Scenario", "subScenario"),
             measure.vars = takeCols)

# Calculate summary stats of variables for each stage across years
temp <- temp %>%
  group_by(Year, Stage, Scenario, subScenario, variable) %>%
  summarise(
    value_mean = mean(value)) %>%
  arrange(variable,Year, Stage) %>%
  droplevels

# Merge all variables
trackParams <- bind_rows(temp, cors) %>% 
  dplyr::mutate(Rep = Rep, .before = Stage)

# Remove objects
rm(temp, cors)
