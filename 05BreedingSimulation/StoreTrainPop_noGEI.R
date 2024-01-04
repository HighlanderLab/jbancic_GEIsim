# Store training population data for GS
#################################################
# Assuming GS starts with 3 years worth of historical records
# The number of records is maintained at 3 years throughout the simulation

#-----------------------------------------------------------------------
# Save MET data from inbred and hybrid trials
#-----------------------------------------------------------------------
if (year == 1) {
  cat(" Start saving training records \n")
  MET = rbind(createMET(pop = PYT, Stage = "PYT",
                        Year = year, nEnvs = 1, nReps = 1),
              createMET(pop = AYT, Stage = "AYT",
                        Year = year, nEnvs = 1, nReps = 1),
              createMET(pop = EYT, Stage = "EYT",
                        Year = year, nEnvs = 1, nReps = 1))
  
} else if (year > 1 && year < (nYearsTP+1)) {
  cat(" Save training records \n")
  # Create initial training population with 4 years of records in burnin
  MET = rbind(MET,
              createMET(pop = PYT, Stage = "PYT",
                        Year = year, nEnvs = 1, nReps = 1),
              createMET(pop = AYT, Stage = "AYT",
                        Year = year, nEnvs = 1, nReps = 1),
              createMET(pop = EYT, Stage = "EYT",
                        Year = year, nEnvs = 1, nReps = 1))
} else if (year > nYearsTP) {
  cat(" Save training records \n")
  # Maintain training population with 4 years of records in future stage
  remove = sum(table(MET$Year)[1])
  MET = rbind(MET[-c(1:remove),],
              createMET(pop = PYT, Stage = "PYT",
                        Year = year, nEnvs = 1, nReps = 1),
              createMET(pop = AYT, Stage = "AYT",
                        Year = year, nEnvs = 1, nReps = 1),
              createMET(pop = EYT, Stage = "EYT",
                        Year = year, nEnvs = 1, nReps = 1))
}

#-----------------------------------------------------------------------
# Save genotypic data of inbreds
#-----------------------------------------------------------------------
if (year == 1) {
  geno_PYT <- pullSnpGeno(PYT)
  geno_AYT <- pullSnpGeno(AYT)
  geno_EYT <- pullSnpGeno(EYT)
  GenoPop  <- rbind(geno_PYT,geno_AYT,geno_EYT)
  rm(geno_PYT,geno_AYT,geno_EYT)
} else if (year > 1) {
  geno_PYT <- pullSnpGeno(PYT)
  geno_AYT <- pullSnpGeno(AYT)
  geno_EYT <- pullSnpGeno(EYT)
  GenoPop  <- rbind(GenoPop,geno_PYT,geno_AYT,geno_EYT)
  GenoPop  <- GenoPop[which(rownames(GenoPop) %in% MET$Id == T),]
  rm(geno_PYT,geno_AYT,geno_EYT)
}