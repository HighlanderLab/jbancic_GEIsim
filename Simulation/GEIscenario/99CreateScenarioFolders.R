## This R script creates separate folder and runs the job for each scenario 
# Launch this script from R on cluster computer

# Import table with all scenarios
testComplete <- read.csv("testComplete.csv")
dim(testComplete)

# Create a folder for each scenario and transfer files
for (i in 1:dim(testComplete)[1]){
  # Create a scenario folder and add the parameters
  dir.create(paste0(i))
  getWD <- getwd()
  setwd(paste0(getWD,"/",i))
  testScenario <- testComplete[i,]
  write.table(testScenario,"testScenario.csv", sep = " ", 
              col.names = F, row.names = F, quote = F)
  # Make copy of files in the scenario folder
  my_files <- c("02RUNMEwithGEI_scenarios.R",
                "AdvanceYear_withGEI.R",
                "CreateCorrelations.R",
                "FillPipeline_withGEI.R",
                "00launchScript.sh") # add launch sh file
  file.copy(from = paste0(getWD,"/", my_files),
            to   = paste0(paste0(getWD,"/",i,"/"), my_files))
  my_files <- c("Functions_asrGEI.R",
                "UpdateParents.R")
  file.copy(from = paste0(getWD,"/../", my_files),
            to   = paste0(paste0(getWD,"/",i,"/"), my_files))
  # Set scenario to run in replicates (call launch sh file) 
  system("./00launchScript.sh")
  # Go back to home directory
  setwd(getWD)
}

# Check which jobs failed and resubmit them
testComplete <- read.csv("testComplete.csv")
reps = 1:10
for (i in 1:dim(testComplete)[1]){
  # Enter scenario folder
  getWD <- getwd()
  setwd(paste0(getWD,"/",i))
  # Get job name 
  testScenario <- read.csv(file = "testScenario.csv",sep = " ",header = F)
  scenarioName = paste0("trackParams_withGEI_",
                        as.integer(testScenario[1]),
                        as.integer(testScenario[2]))
  # Check which jobs are missing 
  temp <- list.files(pattern = scenarioName)
  temp <- gsub(pattern = paste0(scenarioName,"."), replacement = "", temp)
  temp <- gsub(pattern = ".rds", replacement = "", temp)
  temp <- as.numeric(as.character(temp))
  temp <- temp[order(temp)]
  rerun <- reps[!(reps %in% temp)]
  # Rerun missing jobs
  for (i in rerun) {
    system(paste0("qsub -t ", i, " 02RUNMEwithGEI_scenarios.R")) # specify script
  }
  # Go back to home directory
  setwd(getWD)
}

# Collect files from each folder and transfer them to Results folder
dir.create("Results")
for (i in 1:dim(testComplete)[1]){
  # Change directory
  getWD <- getwd()
  setwd(paste0(getWD,"/",i))
  # Move .rds files to Results directory
  listFiles <- list.files(pattern = ".rds")
  file.copy(from = paste0(getwd(),"/", listFiles),
            to = paste0(paste0(getWD,"/Results/"), listFiles))
  # Go back to home directory
  setwd(getWD)
}
