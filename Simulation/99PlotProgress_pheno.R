rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(coda)
library(gridExtra)
library(reshape2)
library(ggpubr)
library(wesanderson)

nReps = 10

# Import table with all scenarios
scenarios <- c("trackParams_noGEI",
               "trackParams_withGEI_0100",
               "trackParams_withGEI_8020",
               "trackParams_withGEI_6040",
               "trackParams_withGEI_4060",
               "trackParams_withGEI_2080",
               "trackParams_withGEI_1000")

# Check if any files missing
temp <- paste(rep(scenarios,each = nReps), ".", rep(1:nReps, times = length(scenarios)), ".rds",sep = "")
temp[!file.exists(temp)]

# Read in output data
rawData = vector("list",nReps*length(scenarios))
i = 0L
for(SCENARIO in scenarios){
  for(REP in c(1:nReps)){
    i = i+1L
    FILE = paste0(SCENARIO,".",REP,".rds")
    temp = readRDS(FILE)
    # temp$Scenario = SCENARIO
    rawData[[i]] = temp
  }
}
rawData = bind_rows(rawData)
colnames(rawData)
str(rawData)

# Specify theme options
optns <- theme_bw(base_size = 16, base_family = "sans") +
  theme(panel.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 14,colour = "black"),
        axis.title = element_text(size = 18),
        axis.title.x=element_blank(),
        axis.text.x=element_text(hjust = .9, angle = 45),
        strip.text = element_text(face = "bold", size = 18, colour = 'black'))

palette <- (c(rep("black",1),
              rep(palette("Set1")[2],1),
              rep("#66CDAA",1),
              rep("#458B74",1),
              rep("#CD8500",1),
              rep("#8B5A00",1),
              rep("#CD6090",1),
              rep("#8B3A62",1)))

####################################################################
# Plot 
####################################################################
temp <- rawData %>% mutate_at(c(1:2), as.numeric)
str(temp)

# Calculate accuracies per stage and melt into one column
cors <- temp %>%
  group_by(Year, Stage, Scenario) %>%
  summarize(acc_obsMean=cor(trueMean, obsMean),
            acc_estMean=cor(trueMean, estMean),
            acc_obsSD=cor(trueSD, obsSD),
            acc_estSD=cor(trueSD, estSD))
cors <- melt(cors,
             id.vars = c("Year","Stage","Scenario"),
             measure.vars = 4:7)
cors

# Melt all other measures into one column
takeCols <- colnames(temp)[5:16]
temp <- melt(temp,
             id.vars = c("Year","Rep","Stage","Scenario"),
             measure.vars = takeCols)

# Calculate summary stats for scenarios across stages and years
temp <- temp %>%
  group_by(Year, Stage, Scenario, variable) %>%
  select(value) %>% 
  summarise(#value = var(value),
    value = mean(value),
    se = sd(value)/sqrt(length(nReps))) %>%
  arrange(variable,Year, Stage) %>%
  filter(variable != "MeanG", variable != "StabilityG") %>%
  droplevels

# Merge all measures
temp <- bind_rows(temp, cors)

# Check data
levels(temp$variable)

# Rename and recode factors
temp$variable <- recode_factor(temp$variable,
                               # `MeanG`    = "MeanG",
                               `trueMean` = "trueMean",
                               `trueSD`   = "trueSD",
                               `trueSI`   = "trueSI",
                               `obsMean`  = "obsMean",
                               `obsSD`    = "obsSD",
                               `obsSI`    = "obsSI",
                               `estMean`  = "estMean",
                               `estSD`    = "estSD",
                               `estSI`    = "estSI",
                               `VarG`     = "VarG",
                               # `StabilityG` = "StabilityG",
                               `acc_obsMean` = "acc_obsMean",
                               `acc_estMean` = "acc_estMean",
                               `acc_obsSD`   = "acc_obsSD",
                               `acc_estSD`   = "acc_estSD")

temp$Stage <- recode_factor(temp$Stage, `HDRW` = "HDRW (estMean)",`PYT` = "PYT (estMean)",`AYT` = "AYT (estSI)", 
                            `EYT` = "EYT (estSI)", `Release` = "Release (estSI)")

temp$Scenario <- recode_factor(temp$Scenario,
                               `noGEI` = "noGEI",
                               `withGEI_1000` = "withGEI_100-0",
                               `withGEI_8020` = "withGEI_80-20",
                               `withGEI_6040` = "withGEI_60-40",
                               `withGEI_4060` = "withGEI_40-60",
                               `withGEI_2080` = "withGEI_20-80",
                               `withGEI_0100` = "withGEI_0-100")

# Plot
(p <- ggplot(temp, aes(x = Year, y = value, col = Scenario)) +
  facet_grid(variable ~ Stage, scales = "free") + 
  geom_line(size = 0.8) +
  # geom_pointrange(aes(ymin = mean-se, ymax = mean+se), size = 1) 
  # # geom_pointrange(temp[temp$colour == "B",], mapping = aes(ymin = mean-se, ymax = mean-se, col = colour)) +
  # # scale_y_continuous("Heterosis", breaks = c(0,10,20, 40, 80)) +
  # scale_y_continuous("Heterosis", limits = c(0,13)) +
  # # scale_x_discrete("Program") +
  scale_colour_manual(values = palette, aesthetics = c("colour", "fill")) +
  optns +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 10)))

ggsave(plot = p, filename ="trackProgress_pheno.pdf",scale = 1.5)
ggsave(plot = p, filename ="trackProgress_pheno.jpg",scale = 1.5)
