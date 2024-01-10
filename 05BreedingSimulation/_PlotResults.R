rm(list=ls())
library(tidyverse)
library(dplyr)
library(coda)
library(gridExtra)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(wesanderson)

# Specify theme options
optns <- theme_bw(base_size = 16, base_family = "sans") +
  theme(panel.background = element_blank(),
        # legend.title = element_blank(),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text  = element_text(size = 14),
        legend.position = "right",
        plot.title = element_text(size = 15, face = "bold"),
        axis.text  = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 18),
        # axis.title.x = element_blank(),
        axis.text.x  = element_text(hjust = .9, angle = 45),
        strip.text   = element_text(face = "bold", size = 18, colour = 'black'))

####################################################################
# Process dataset 
####################################################################
nReps <- 1

temp <- list.files(pattern = "Results.*GEI") %>%
  # temp <- list.files(pattern = "Low.*rds") %>%
  map(readRDS) %>% 
  bind_rows()

# Check if all is okay
table(temp$Scenario,temp$subScenario)
# table(temp$Rep, temp$Scenario, temp$subScenario)

temp$Stage <- as.factor(temp$Stage)
# temp$Scenario <- as.factor(temp$Scenario)
# temp$subScenario <- as.factor(temp$subScenario)
temp$variable <- as.factor(temp$variable)
# Remove first 3 years because they serve as burnin
temp$Year <- temp$Year - 3

# Calculate summary stats for scenarios across stages and years
temp <- temp %>%
  group_by(Year, Stage, Scenario, subScenario, variable) %>%
  summarise(
    value_mean = mean(value_mean),
    value_se   = sd(value_mean, na.rm = T)/sqrt(nReps)) %>%
  arrange(variable,Year, Stage) %>%
  # filter(variable != "MeanG") %>%
  droplevels

# Set MET to TPE accuracy of noGEI scenario to 1
temp[temp$variable %in% c("acc_Mean_METtoTPE") & temp$Scenario %in% c("Pheno_noGEI","GS_noGEI"),]$value_mean <- 1
temp[temp$variable %in% c("acc_Mean_subMETtoTPE") & temp$Scenario %in% c("Pheno_noGEI","GS_noGEI"),]$value_mean <- 1
# Add noGEI scenario to TPE as well
temp[temp$variable %in% c("Mean_TPE") & temp$Scenario %in% c("Pheno_noGEI","GS_noGEI"),]$value_mean <-
  temp[temp$variable %in% c("Mean_MET") & temp$Scenario %in% c("Pheno_noGEI","GS_noGEI"),]$value_mean
temp[temp$variable %in% c("VarG_TPE") & temp$Scenario %in% c("Pheno_noGEI","GS_noGEI"),]$value_mean <-
  temp[temp$variable %in% c("VarG_MET") & temp$Scenario %in% c("Pheno_noGEI","GS_noGEI"),]$value_mean
temp[temp$variable %in% c("acc_Mean_estMETtoTPE") & temp$Scenario %in% c("Pheno_noGEI","GS_noGEI"),]$value_mean <-
  temp[temp$variable %in% c("acc_Mean_estMETtoMET") & temp$Scenario %in% c("Pheno_noGEI","GS_noGEI"),]$value_mean
# Set mean and variance for Pheno_GEI programs because MET is limited to single year
temp[temp$variable %in% c("Mean_MET") & temp$Scenario %in% c("PhenoGEI_1-0"),]$value_mean <-
  temp[temp$variable %in% c("estMean_MET") & temp$Scenario %in% c("PhenoGEI_1-0"),]$value_mean
temp[temp$variable %in% c("VarG_MET") & temp$Scenario %in% c("PhenoGEI_1-0"),]$value_mean <-
  temp[temp$variable %in% c("VarG_subMET") & temp$Scenario %in% c("PhenoGEI_1-0"),]$value_mean
# Set MET-TPE alignment for Pheno_GEI programs because MET is limited to single year
temp[temp$variable %in% c("acc_Mean_METtoTPE") & temp$Scenario %in% c("PhenoGEI_1-0"),]$value_mean <-
  temp[temp$variable %in% c("acc_Mean_subMETtoTPE") & temp$Scenario %in% c("PhenoGEI_1-0"),]$value_mean
# Set MET-TPE alignment for EYT stage where GS is not used
temp[temp$variable %in% c("acc_Mean_METtoTPE") & temp$Stage %in% c("EYT"),]$value_mean <-
  temp[temp$variable %in% c("acc_Mean_subMETtoTPE") & temp$Stage %in% c("EYT"),]$value_mean

# Rename and recode variables
temp$variable <- recode_factor(temp$variable,
                               `Mean_TPE` = "Main effects in TPE",
                               `SD_TPE`   = "SD_TPE",
                               `SI_TPE`   = "SI_TPE",
                               `Mean_MET` = "Main effects in MET",
                               `SD_MET`   = "SD_MET",
                               `SI_MET`   = "SI_MET",
                               `Mean_subMET`  = "Mean_subMET",
                               `SD_subMET`    = "SD_subMET",
                               `SI_subMET`    = "SI_subMET",
                               `estMean_MET`  = "estMean_MET",
                               `estSD_MET`    = "estSD_MET",
                               `estSI_MET`    = "estSI_MET",
                               `MeanG_TPE`    = "MeanG_TPE",
                               `MeanG_subTPE` = "MeanG_subTPE",
                               `VarG_TPE`     = "Main effects in TPE ",
                               `VarG_subTPE`  = "VarG_subTPE",
                               `VarMeanG_TPE` = "VarMeanG_TPE",
                               `VarMeanG_subTPE` = "VarMeanG_subTPE",
                               `MeanA_TPE`    = "MeanA_TPE",
                               `MeanA_subTPE` = "MeanA_subTPE",
                               `VarA_TPE`     = "VarA_TPE",
                               `VarA_subTPE`  = "VarA_subTPE",
                               `MeanD_TPE`    = "MeanD_TPE",
                               `MeanD_subTPE` = "MeanD_subTPE",
                               `VarD_TPE`     = "VarD_TPE",
                               `VarD_subTPE`  = "VarD_subTPE",
                               `MeanG_MET`    = "MeanG_MET",
                               `MeanG_subMET` = "MeanG_subMET",
                               `VarG_MET`     = "Main effects in MET ",
                               `VarG_subMET`  = "VarG_subMET",
                               `VarMeanG_MET` = "VarMeanG_MET",
                               `VarMeanG_subMET` = "VarMeanG_subMET",
                               `MeanA_MET`    = "MeanA_MET",
                               `MeanA_subMET` = "MeanA_subMET",
                               `VarA_MET`     = "VarA_MET",
                               `VarA_subMET`  = "VarA_subMET",
                               `MeanD_MET`    = "MeanD_MET",
                               `MeanD_subMET` = "MeanD_subMET",
                               `VarD_MET`     = "VarD_MET",
                               `VarD_subMET`  = "VarD_subMET",
                               `estVar_MET`   = "estVar_MET",
                               `acc_Mean_subMETtoMET` = "Mean_subMETtoMET", # cor (gv,gv)
                               `acc_Mean_estMETtoTPE` = " Main effects in TPE",
                               `acc_Mean_estMETtoMET` = " Main effects in MET",
                               `acc_Mean_METtoTPE`    = " MET-TPE alignment", # cor (gv,gv)
                               `acc_Mean_subMETtoTPE` = "Mean_subMETtoTPE", # cor (gv,gv)
                               `acc_SD_METtoTPE`      = "SD_METtoTPE",
                               `acc_SD_subMETtoTPE`   = "SD_subMETtoTPE",
                               `acc_SD_subMETtoMET`   = "SD_subMETtoMET",
                               `acc_SD_estMETtoMET`   = "SD_estMETtoMET",
                               `acc_SD_setMETtoTPE`   = "SD_setMETtoTPE")
levels(temp$variable)

# Rename and recode stages
temp$Stage <- recode_factor(temp$Stage,
                            `Parents` = "Parents",
                            `HDRW`    = "HDRW",
                            `postHDRW`    = "postHDRW",
                            `PYT`     = "PYT",
                            `AYT`     = "AYT",
                            `EYT`     = "EYT",
                            `Release` = "Release")
levels(factor(temp$Stage))

# Rename and recode subScenarios - i.e. GEI complexity
temp$GEI <- temp$subScenario 
temp$GEI[temp$Scenario == "Pheno_noGEI"] <- "No"
temp$GEI[temp$Scenario == "GS_noGEI"] <- "No"
temp$GEI <- recode_factor(temp$GEI,
                          `no GEI`      = "No",
                          `LowGEI`      = "Low",
                          `ModerateGEI` = "Moderate",
                          `HighGEI`     = "High")
levels(factor(temp$GEI))

# Create factor for Program for investigation
levels(factor(temp$Scenario))
temp <- temp %>%
  filter((Scenario %in% c("Pheno_noGEI",
                          "GS_noGEI",
                          "PhenoGEI_1-0",
                          "GSGEI_1-0"
                          ))) %>% droplevels() 

temp$Program <- recode_factor(temp$Scenario,
                              `Pheno_noGEI`  = "Phenotypic Selection",
                              `GS_noGEI`     = "Genomic Selection",
                              `PhenoGEI_1-0` = "Phenotypic Selection",
                              `GSGEI_1-0`    = "Genomic Selection"
)
levels(factor(temp$Program))

# Colour of lines
palette <- (c(rep("black",1),
              rep("#66CDAA",1),
              rep("#458B74",1),
              rep("#8B4513",1),
              rep(palette("Set1")[2],1)
))

#-----------------------------------------------------------------------
# Plot - mean values
#-----------------------------------------------------------------------
(a <- temp %>%
    filter((variable %in% c("Main effects in TPE","Main effects in MET"))) %>%
    filter((Stage %in% c("HDRW"))) %>%
    # filter(Year > -3 & Year < 18) %>%
    droplevels() %>% 
    ggplot(aes(x = Year, y = value_mean)) +
    facet_grid(variable ~ Program, scales = "fixed") +
    geom_line(aes(color = GEI), linewidth = 0.8) +
    # facet_grid(Stage ~ variable, scales = "fixed") + # plot all stages
    # geom_line(aes(linetype = GEI, color = Scenario), linewidth = 0.8) +  # plot all stages
    scale_linetype_manual(values=c("twodash", "solid")) +
    ylab("Genetic gain") + 
    scale_y_continuous(limits = c(-0.5, 12), breaks = seq(0, 12, 2)) +
    xlim(0,20) +
    scale_colour_manual(values = palette, aesthetics = c("colour", "fill")) +
    optns +
    theme(strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          # strip.text.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm")))
ggsave(plot = a, filename ="01trackProgress_Mean.pdf", width = 4, height = 3, scale = 2)

##
(b <- temp %>%
    filter((variable %in% c("Main effects in TPE ", "Main effects in MET "))) %>%
    filter((Stage %in% c("HDRW"))) %>%
    # filter(Year > 20) %>%
    droplevels() %>% 
    ggplot(aes(x = Year, y = value_mean)) +
    facet_grid(variable ~ Program, scales = "fixed") +
    geom_line(aes(color = GEI), linewidth = 0.8) +
    # facet_grid(Stage ~ variable, scales = "fixed") + # plot all stages
    # geom_line(aes(linetype = GEI, color = Scenario), linewidth = 0.8) +  # plot all stages
    scale_linetype_manual(values=c("twodash", "solid")) +
    ylab("Variance") + 
    xlim(0,20) +
    scale_colour_manual(values = palette, aesthetics = c("colour", "fill")) +
    optns +
    theme(strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          # strip.text.y = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm")))
ggsave(plot = b, filename ="02trackProgress_Var.pdf", width = 4, height = 3, scale = 2)

##
(c <- temp %>%
    filter((variable %in% c(" MET-TPE alignment"," Main effects in TPE"," Main effects in MET"))) %>%
    filter((Stage %in% c("HDRW"))) %>%
    # filter(Year > 20) %>%
    droplevels() %>% 
    ggplot(aes(x = Year, y = value_mean)) + 
    facet_grid(variable ~ Program, scales = "fixed") +
    geom_line(aes(color = GEI), linewidth = 0.8) +
    # facet_grid(Stage ~ variable, scales = "fixed") + # plot all stages
    # geom_line(aes(linetype = GEI, color = Scenario), linewidth = 0.8) +  # plot all stages
    scale_linetype_manual(values=c("twodash", "solid")) +
    ylab("Accuracy") +
    ylim(0,1) +
    xlim(0,20) +
    scale_colour_manual(values = palette, aesthetics = c("colour", "fill")) +
    optns +
    theme(strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          plot.margin = unit(c(0, 0, 0, 0), "cm")))
ggsave(plot = c, filename ="03trackProgress_Acc.pdf", width = 5, height = 5, scale = 1.7)
