# Script removes 10 oldest parents each year and replaces them 
# with most recent genotypes from the EYT stage to form a new crossing block

# Drop 10 parents
Parents = Parents[11:nParents]

# Update with new 10 parents from the EYT stage
Parents = c(Parents, EYT)

